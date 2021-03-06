﻿/*---------------------------------------------------------------*
 *   Name: qsim.sas                                        	     *
 *   Title: Simulation of queueing systems						 *
 *---------------------------------------------------------------*
 * Author:  Dmitry Shopin            <dmitry.shopin@gmail.com>   *
 * Created: 7 Sep 2017      		                             *
 * 							                                     *
 * Version: 1.0                                                  *
 *---------------------------------------------------------------*/
 /*
Description:

 The QSIM macro carries out simulation of a queueing system  with specified
number of servers and distirbutions of interarrival and service times.
The assumptions are:
	-infinite population size
	-arrival process independent of the state of the system
	-infinite system capacity (infinite queue)

(Or G/G/s/FIFO/Inf/Inf in Kendall notation)

Usage:

 The QSIM macro is called with keyword parameters.
 The arguments may be listed within parentheses in any order, separated by commas. For example:

   %qsim(ntask=10000, nserv=6, iatdist=exponential(0.04), servdist=exponential(0.007), seed=-1)


==Parameters:


* NTASK=	Specifies the number of simulated tasks. Default=1000.

* NSERV=	Specifies the number of servers. Default=1

* IATDIST=  Specifies the name of the distribution of interarrival time and corresponding parameters
			in the form: DISTRIBUTION_NAME(par1, par2, ...). Parameters must be specified as following:
			BERNOULLI(p)
			BETA(a,b)
			BINOMIAL(p,n)
			CHISQUARE(df)
			ERLANG(shape, scale) (scale=1/rate)
			EXPONENTIAL(lambda)
			F(n,d)
			GAMMA(shape, scale)
			GEOMETRIC(p)
			HYPERGEOMETRIC(N,R,n)
			LOGNORMAL(mu, sigma)
			NEGBINOMIAL(p,k)
			PARETO(alpha, sigma)
			POISSON(lambda)
			TABLE(p1,p2,...,pn)
			TRIANGULAR(min, max, mode)
			UNIFORM(min, max)
			WEIBULL(shape, scale)

			No default.

* SERVDIST=  Specifies the name of the distribution of service time and corresponding parameters same as IATDIST

* SEED=		Seed value for random number generations. Default=-1 (generating seeds from the system clock)

 =

Output:

*ODS report with some calculated metrics and plots

*Dataset SIMTASKS, where each row contains data for one task (arrival, start time of serving etc)

*Dataset SIMSERVERS with events for each individual server ('engage'=start serving, 'release' = getting available again)

*Dataset SIMQUEUE storing the queue length at any point in time (each row is a period of constant queue length)

*Dataset OUTSTAT with some calculated metrics in one row
*/


%macro qsim(ntask=1000,nserv=1,iatdist=,servdist=,seed=-1);


	/*Macro to output error message if IATDIST or SERVDIST are invalid*/
	%macro distinvalmsg;
		%put ERROR: %upcase(&key) argument must have form Distribution_Name(param-1, ..., param-n);
		%put ERROR: %upcase(&key) = &&&key;
		%put NOTE: Valid distribution names are BERNOULLI, BETA, BINOMIAL, CHISQUARE, ERLANG, EXPONENTIAL, F, GAMMA, GEOMETRIC, HYPERGEOMETRIC, LOGNORMAL, NEGBINOMIAL, PARETO, POISSON, TABLE, TRIANGULAR, UNIFORM, WEIBULL;
		%abort;
	%mend distinvalmsg;


	/*Macro to generate various random variables*/
	%macro randgen(key=, varname=);
		%if &&&key.name=PARETO %then %do;
			U = rand("Uniform"); /*Pareto is not supported by RAND(), so we use inverse transform sampling*/
			&varname = -  &&&key.parm2/ &&&key.parm1 * (U** &&&key.parm1 - 1);
		%end;

		/*Some functions are supported by RAND() in standard form only. Shift and scale them*/
		%else %if &&&key.name=EXPONENTIAL %then
			&varname = rand("EXPONENTIAL") / &&&key.parm1; /*If x~Exp(1) then x/lambda ~ Exp(lambda)*/

		%else %if &&&key.name=LOGNORMAL %then /*if x~N(mu,var) then exp(x)~Lognorm(mu,var)*/
			&varname = exp(rand("NORMAL",  &&&key.parm1,  &&&key.parm2));

		%else %if &&&key.name = GAMMA %then /*if x~Gamma(shape,1) then scale*x~Gamma(shape, scale) */
			&varname = rand("GAMMA",  &&&key.parm1) *  &&&key.parm2;

		%else %if &&&key.name = ERLANG %then /*if x~Erlang(shape,1) then scale*x~Gamma(shape, scale) */
			&varname = rand("ERLANG",  &&&key.parm1) *  &&&key.parm2;

		%else %if &&&key.name=TRIANGLE %then /*if x~Tri(0,1,(mode-min)/(max-min)) then a+(b-a)*x~Tri(min, max, mode) */
			&varname = &key.parm1 + ( &&&key.parm2- &&&key.parm1)*rand("TRIANGLE",( &&&key.parm3- &&&key.parm1)/( &&&key.parm2- &&&key.parm1));

		%else %if &&&key.name=UNIFORM %then /*if x~U(0,1) then a+(b-a)*x~U(min, max) */
			&varname =  &&&key.parm1 + ( &&&key.parm2- &&&key.parm1)*rand("UNIFORM");

		%else &varname = rand("&&&key.name", &&&key.parms);;
	%mend randgen;



	/*Macro for pre-processing IATDIST and SERVDIST keywords*/
	%macro distproc(key) / minoperator;
		/*Extract distribution and parameters for IAT*/
		%global &key.name &key.parms;
		%let firstbracket = %index(&&&key,%str(%());
		%let lastbracket = %index(&&&key,%str(%)));

		/*Check the structure of IATDIST and SERVDIST*/
		%if &firstbracket < 5 or %eval(&lastbracket - &firstbracket)<2 %then %distinvalmsg;

		%let &key.name = %upcase(%substr(&&&key,1,%eval(&firstbracket-1)));
		%let &key.parms=%substr(&&&key, %eval(&firstbracket+1),%eval(&lastbracket-&firstbracket-1));

		/*Check the names of IATDIST and SERVDIST*/
		%if not(%substr(&&&key.name,1,4) in BERN BETA BINO CHIS ERLA EXPO F GAMM GEOM HYPE LOGN NEGB POIS TABL TRIA UNIF WEIB PARE)
			%then %distinvalmsg;

		/*unify names*/
		%if %substr(&&&key.name,1,4)=BERN %then %let &key.name=BERNOULLI;
		%else %if %substr(&&&key.name,1,4)=BINO %then %let &key.name=BINOMIAL;
		%else %if %substr(&&&key.name,1,4)=CHIS %then %let &key.name=CHISQUARE;
		%else %if %substr(&&&key.name,1,4)=ERLA %then %let &key.name=ERLANG;
		%else %if %substr(&&&key.name,1,4)=EXPO %then %let &key.name=EXPONENTIAL;
		%else %if %substr(&&&key.name,1,4)=GAMM %then %let &key.name=GAMMA;
		%else %if %substr(&&&key.name,1,4)=GEOM %then %let &key.name=GEOMETRIC;
		%else %if %substr(&&&key.name,1,4)=HYPE %then %let &key.name=HYPERGEOMETRIC;
		%else %if %substr(&&&key.name,1,4)=LOGN %then %let &key.name=LOGNORMAL;
		%else %if %substr(&&&key.name,1,4)=NEGB %then %let &key.name=NEGBINOMIAL;
		%else %if %substr(&&&key.name,1,4)=POIS %then %let &key.name=POISSON;
		%else %if %substr(&&&key.name,1,4)=TABL %then %let &key.name=TABLE;
		%else %if %substr(&&&key.name,1,4)=TRIA %then %let &key.name=TRIANGLE;
		%else %if %substr(&&&key.name,1,4)=UNIF %then %let &key.name=UNIFORM;
		%else %if %substr(&&&key.name,1,4)=WEIB %then %let &key.name=WEIBULL;
		%else %if %substr(&&&key.name,1,4)=PARE %then %let &key.name=PARETO;

		/*Extract individual parameters*/
		%do i=1 %to %sysfunc(countw("&&&key.parms",%str(,)));
			%global &key.parm&i;
			%let &key.parm&i=%scan(%bquote(&&&key.parms),&i,%str(,));
		%end;

		/*Checks if number of parameters is valid. Extra parameters are ignored*/
		%if %sysfunc(countw(%bquote(&&&key.parms),%str(,)))<3 and  &&&key.name in (HYPERGEOMETRIC TRIANGLE)
			or %sysfunc(countw(%bquote(&&&key.parms),%str(,)))<2 and &&&key.name in (BETA BINOMIAL ERLANG F GAMMA LOGNORMAL NEGBINOMIAL UNIFORM WEIBULL)
			or %sysfunc(countw(%bquote(&&&key.parms),%str(,)))<1
		%then %do;
			%put ERROR: Not enough parameters for &&&key.name distribution; %abort;
		%end;

		/*Check if parameters are valid for some distributions by testing the function*/
		data _null_;
			%randgen(key=&key, varname=rnum);
			call symputx("rnum",rnum);
		run;
		%if &rnum=. %then %do;
			%put ERROR: Invalid parameters for &&&key.name distribution; %abort;
		%end;

	%mend distproc;
	%distproc(IATDIST)
	%distproc(SERVDIST)



	/*Main simulation process*/
	data 	simqueue(keep=event clock)
			simservers(keep=server_id event clock task_id)
			simtasks(keep=task_id event clock);

		length enqueue_time service_time task_id idle_time server_id release_time 8;
		length event $8.;

		call streaminit(&seed);
		
		if _N_=1 then do;	
			/*hash object for the queue*/
			dcl hash Q(ordered:'a');
			dcl hiter iQ('Q');
			Q.defineKey('enqueue_time');
			Q.defineData('enqueue_time', 'service_time', 'task_id');
			Q.defineDone();

			/*hash object for the current state of each server*/
			dcl hash S(ordered:'a');
			dcl hiter iS('S');
			S.defineKey('server_id');
			S.defineData('server_id', 'idle_time', 'release_time', 'task_id');
			S.defineDone();

			/*create neccessary number of servers*/
			idle_time = 0;
			do server_id=1 to &nserv;
				rc0=S.add();
			end;
		end;
		
		clock = 0;
		
		do next_task_id = 1 to &ntask; /*Main loop, each iteration - 1 new task*/

			/*Generate interarrival time for the next call*/
			%randgen(key=IATDIST, varname=IAT)
			next_task = clock + IAT;

			/*Generate service time for the next call*/
			%randgen(key=SERVDIST, varname=next_service_time)



			do until(clock=next_task); /*start of loop for events before arrival of the next task*/

				 /********************************************************/
				/*Transferring tasks from the queue to available servers*/
			   /********************************************************/

				do while(Q.num_items > 0); /*only if there is at least one task in the queue*/
					rc = iS.first();/*check if there are available servers and who's idle the longest*/
					maxidle = .;
					maxidleid = .;
					do while(rc = 0);
						if (idle_time^=. and maxidle^=. and idle_time > maxidle) or (idle_time^=. and maxidle=.) then do;
							maxidle = idle_time;
							maxidleid = server_id;
						end;
						rc = iS.next();
					end;
					server_id = maxidleid;				
					rc = S.find();
					if rc = 0 then do;
						/*dequeue task*/
						iQ.first();
						event = 'dequeue';
						output simqueue;
						rc0 = iQ.delete(); /*have to delete iterator before removing element*/
						rc0 = Q.remove();
						iQ = _new_ hiter('Q'); /*recreating iterator*/

						/*engage server*/
						event = 'engage';
						output simservers;	
						event = 'start';
						output simtasks;
						idle_time = .;
						release_time = clock + service_time;
						rc0 = S.replace();
					end;
					else leave;
				end;


				 /************************************************************/
				/*Search for the closest server release before the next task*/
			   /************************************************************/
				rc = iS.first();
				nextrelease = .;
				nextreleaseid = .;
				do while(rc = 0);
					if (release_time^=. and nextrelease^=. and release_time < nextrelease) or (release_time^=. and nextrelease=.) then do;
						nextrelease = release_time;
						nextreleaseid = server_id;
					end;
					rc = iS.next();
				end;
				server_id = nextreleaseid;				
				rc = S.find();
				if rc = 0 and release_time < next_task then do;
					jump = release_time - clock; /*we need to store this jump to update servers idle time below*/ 
					clock = release_time; /*jumping to the next release time*/
				end;
				else do; /*if no releases before the next task then jump to the next task*/ 
					jump = next_task - clock;
					clock = next_task; 
				end;

				/*releasing server(s)*/
				rc = iS.first();
				do while(rc = 0);
					if release_time = clock then do;
						event = 'release';
						output simservers;
						event = 'end';
						output simtasks;
						/*update server hash record*/
						idle_time = 0;
						release_time = .;
						task_id = .;
						rc0=S.replace();
					end;
					rc=iS.next();
				end;

				/*update idle times for servers*/
				rc = iS.first();
				do while(rc = 0);
					if release_time = . then do;
						idle_time + jump;
						rc0 = S.replace();
					end;
					rc = iS.next();
				end;
			end; /*end of loop for events before arrival of the next task*/

			/*processing new call*/
			task_id = next_task_id;
			event = 'enqueue';
			output simqueue;

			event = 'arrival';
			output simtasks;
			enqueue_time = clock;
			service_time = next_service_time;
			rc=Q.add();
		end;

	run;

	 /***************/	
	/*Tasks dataset*/
   /***************/
	proc sort data=simtasks; by task_id; run;
	proc transpose data=simtasks out=simtasks(drop=_:);
		by task_id;
		id event;
		var clock;
	run;
	data simtasks;
		set simtasks;
		waiting_time = start - arrival;
		service_time = end - start;
		system_time = end - arrival;
		iat = coalesce(arrival - lag(arrival), arrival);
	run;


	   /***************/	
	  /*Queue dataset*/
	 /***************/
	data simqueue;
		set simqueue;
		change = (event='enqueue') - (event='dequeue');
	run;

	data simqueue;
		set simqueue;
		by clock;
		if FIRST.clock then cumchange=0;
		cumchange+change;
		if LAST.clock then output;
		keep clock cumchange;
		rename cumchange=change;
	run;



	/*excl. 0 changes (appeared because of tasks enqueued and dequeued same moment*/
	data simqueue;
		set simqueue;
		if _N_=1 or change^=0;
	run;

	/*count current queue length for each row and time to the next change*/
	data simqueue;
		set simqueue;
		set simqueue(firstobs=2 keep=clock rename=(clock=clock_next));
		L + change;
		duration = clock_next - clock;
		drop clock_next change;
		label L="Queue length"
			  duration="Period of constant queue";
	run;


	   /*****************/	
	  /*Servers dataset*/
	 /*****************/
	proc sort data=simservers; by server_id task_id; run;

	data simservers;
		set simservers;
		retain start end;
		by server_id task_id;
		if event='engage' then start=clock;
		else end=clock;
		duration=end-start;
		if LAST.task_id then output;
		drop event clock;
	run;





	   /***************/	
	  /*Various stats*/
	 /***************/
	proc sql noprint;
		select avg(iat) into :iat_mean from simtasks;
		select avg(service_time) into :service_mean from simtasks;
		select avg(waiting_time) into :wait_mean from simtasks;
		select sum(L*duration)/sum(duration) into :L_mean from simqueue;
		select 1-sum(duration)/(max(start)*&nserv) format=percent8.3 into :Idle from simservers 
	quit;
		

	data _null_;
		file print;
		title 'Simulation Results';
		put / @10 "Interarrival time distribution: &iatdistname (&iatdistparms)"
		    / @10 "Service time distribution: &servdistname (&servdistparms)"
			/ @10 "Number of servers: &nserv"
			/ @10 "Number of tasks: &ntask"
			/ @10 "Seed: &seed"
			// @10 "Interarrival time, mean: &iat_mean"
			/ @10 "Interarrival rate(lambda): %sysevalf(1/&iat_mean)"
			/ @10 "Service time, mean: &service_mean"
			/ @10 "Service rate (mu): %sysevalf(1/&service_mean)"
			/ @10 "Traffic intensity (ro): %sysevalf(&service_mean/(&iat_mean*&nserv))"
			/ @10 "Queue length, mean (Lq): &L_mean"
			/ @10 "Calculated queue length (Little's Law: Lq=lambda*Wq): %sysevalf(%sysevalf(1/&iat_mean)*&wait_mean)"
			/ @10 "Waiting time, mean (Wq): &wait_mean"
			/ @10 "Servers' idle time percent, mean: &Idle"
	;
	run;
	title;

	/*Creating 1-row dataset with the same stats*/
	data outstat;
		nserv=&nserv;
		ntask=&ntask;
		lambda=%sysevalf(1/&iat_mean);
		service_mean=&service_mean;
		mu=%sysevalf(1/&service_mean);
		ro=%sysevalf(&service_mean/(&iat_mean*&nserv));
		IAT=&iat_mean;
		Lq=&L_mean;
		Wq=&wait_mean;
		label
			nserv="Number of servers"
			ntask="Number of tasks"
			IAT="Interarrival time, mean"
			lambda="Interarrival rate"
			service_mean="Service time, mean"
			mu="Service rate"
			ro="Traffic intensity"
			Lq="Queue length, mean"
			Wq="Waiting time, mean";
	run;

	   /***************/	
	  /*    Plots    */
	 /***************/

	proc univariate data=simtasks noprint;
		var waiting_time;
		output out=pctl  pctlpts=50 80 90 95 99 pctlpre=P;
	run;
	data _null_;
		set pctl;
		call symputx("P50", P50);
		call symputx("P80", P80);
		call symputx("P90", P90);
		call symputx("P95", P95);
		call symputx("P99", P99);
	run;

	proc sgplot data=simtasks;
		title 'Interarrival Time Distribution';
		histogram iat;
	run;

	proc sgplot data=simtasks;
		title 'Service Time Distribution';
		histogram service_time;
	run;

	proc sgplot data=simqueue;
		title 'Instantaneous Queue';
		series x=clock y=L;
	run;
	title;

	proc sgplot data=simtasks;
		title 'Waiting Time Distribution';
		histogram waiting_time;
		inset "Mean: &wait_mean"
				"Median: &P50"
				"80% below &P80"
				"90% below &P90"
				"95% below &P95"
				"99% below &P99"
				 / position=topright valuealign=left;
	run;
	title;
	proc datasets nolist; delete pctl; run;

%mend qsim;