/*-------------------------------------------------------------------*
           *    Name: qsim.sas                                        	      *
           *   Title: Simulation of queueing systems							  *
           *-------------------------------------------------------------------*
           *  Author:  Dmitry Shopin            <dmitry.shopin@gmail.com>      *
           * Created: 16 Aug 2017      		                                  *
  * 							                                      *
  * Version: 1.0                                                      *
  *-------------------------------------------------------------------*/
 /*
Description:

 The QSIM macro carries out simulations of queueing systems with user-defined
number of servers and distirbutions of interarrival and service times.
The assumptions are:
	-infinite population size
	-arrival process independent of the state of the system
	-infinite system capacity (infinite queue)

Usage:

 The QSIM macro is called with keyword parameters.
 The arguments may be listed within parentheses in any order, separated by commas. For example:

   %qsim(ntask=10000, nserv=6, iatdist=exponential(0.04), servdist=exponential(0.007), seed=-1)


==Parameters:


* NTASK=	Specifies the number of simulated tasks (jobs, calls etc). Default=1000.

* NSERV=	Specifies the number of servers. Default=1


* IATDIST=  Specifies the name of the distribution of interarrival time and corresponding parameters
			in the form: DISTRIBUTION_NAME(par1, par2, ...). Parameters must be specified in an order
			defined for RAND() function
			(see http://support.sas.com/documentation/cdl/en/lefunctionsref/69762/HTML/default/viewer.htm#p0fpeei0opypg8n1b06qe4r040lv.htm)
			No default.

* SERVDIST=  Same as for DISTARR

* SEED=		Seed value for random number generations. Default=-1 (generation seeds from the system clock)

 =*/


%macro qsim(
  ntask=1000
 ,nserv=1
 ,iatdist=
 ,servdist=
 ,seed=-1
);

options mprint;
/*extract distribution and parameters for IAT*/
%let iatdistname=%substr(%upcase(&iatdist),1,%eval(%index(&iatdist,%str(%())-1));
%let firstbrac = %index(&iatdist,%str(%());
%let lastbrac = %index(&iatdist,%str(%)));
%let iatparms=%substr(%upcase(&iatdist), %eval(&firstbrac+1),%eval(&lastbrac-&firstbrac-1));

%do i=1 %to %sysfunc(countw("&iatparms",%str(,)));
	%let iatparm&i=%scan(%nrbquote(&iatparms),&i,%str(,));
%end;


/*extract distribution and parameters for Service time*/
%let servdistname=%substr(%upcase(&servdist),1,%eval(%index(&servdist,%str(%())-1));
%let firstbrac = %index(&servdist,%str(%());
%let lastbrac = %index(&servdist,%str(%)));
%let servparms=%substr(%upcase(&servdist), %eval(&firstbrac+1),%eval(&lastbrac-&firstbrac-1));

%do i=1 %to %sysfunc(countw("&servparms",%str(,)));
	%let servparm&i=%scan(%nrbquote(&servparms),&i,%str(,));
%end;


data 	simqueue(keep=event clock)
		simservers(keep=server_id event clock)
		simtasks(keep=task_id event clock);

	length enqueue_time service_time task_id idle_time server_id release_time 8;
	length event $8.;

	call streaminit(&seed);
	
	if _N_=1 then do;	
		dcl hash Q(ordered:'a');
		dcl hiter iQ('Q');
		Q.defineKey('enqueue_time');
		Q.defineData('enqueue_time', 'service_time', 'task_id');
		Q.defineDone();

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
	
	do next_task_id = 1 to &ntask;

		/*when next call*/
	%if &iatdistname=PARETO %then %do;
		U = rand("Uniform");
		IAT = - &iatparm2/&iatparm1 * (U**&iatparm1 - 1);
		IAT = ceil(IAT);
		next_task = clock + IAT;
	%end;
	%else %if &iatdistname=EXPONENTIAL %then %do;
		IAT = rand("&iatdistname") / &iatparm1;
		IAT = ceil(IAT);
		next_task = clock + IAT;
	%end;
	%else %do;
		IAT = rand("&iatdistname", &iatparms);
		IAT = ceil(IAT);
		next_task = clock + IAT;
	%end;


		/*generate service time*/
	%if &servdistname=PARETO %then %do;
		U = rand("Uniform");
		next_service_time = - &servparm2/&servparm1 * (U**&servparm1 - 1);
		next_service_time = ceil(next_service_time);
	%end;
	%else %if &servdistname=EXPONENTIAL %then %do;
		next_service_time = rand("&servdistname") / &servparm1;
		next_service_time = ceil(next_service_time);
	%end;
	%else %do;
		next_service_time = rand("&servdistname", &servparms);
		next_service_time = ceil(next_service_time);
	%end;

/*		put 'Call #' next_task_id 'will arrive at ' next_task 'and its duration ' next_service_time;*/

		do clock=clock to next_task-1;
			/*releasing servers*/
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
/*					put 'Call #' task_id 'finished by server #' server_id 'at ' clock; */
				end;
				rc=iS.next();
			end;

			/*dequeueing*/
			do while(Q.num_items > 0);

				/*check if there are available servers and who's idlest*/
				rc = iS.first();
				idlest = .;
				maxidle = .;
				do while(rc = 0);
					if idle_time > maxidle then do;
						maxidle = idle_time;
						idlest = server_id;
					end;
					rc = iS.next();
				end;
				server_id = idlest;
				
				rc = S.find();
				if rc = 0 then do;
					/*dequeue*/
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

/*					put 'Call #' task_id 'dequeued and transferred to server #' server_id 'at ' clock;*/
				end;
				else leave;
			end;

			/*update idle times for servers*/
			rc = iS.first();
			do while(rc = 0);
				if release_time = . then do;
					idle_time + 1;
					rc0 = S.replace();
				end;
				rc = iS.next();
			end;
		end;

		/*processing new call*/
		task_id = next_task_id;
		event = 'enqueue';
		output simqueue;
/*		put 'Call #' task_id ' inqueued at ' clock;*/

		event = 'arrival';
		output simtasks;
		enqueue_time = clock;
		service_time = next_service_time;
		rc=Q.add();
	end;

run;


/*Call stats*/
proc sort data=simtasks; by task_id; run;
proc transpose data=simtasks out=simtasks;
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



/*Queue stats*/
data simqueue;
	set simqueue;
	change = (event='enqueue') - (event='dequeue');
run;
proc freq data=simqueue noprint;	
	tables clock / out=simqueue(drop=percent);
	weight change;
run;
data simqueue;
	set simqueue;
	set simqueue(firstobs=2 keep=clock rename=(clock=clock_next));
	L + count;
	duration = clock_next - clock;
run;

proc sql noprint;
	select avg(iat) into :iat_mean from simtasks;
	select avg(service_time) into :service_mean from simtasks;
	select avg(waiting_time) into :wait_mean from simtasks;
	select sum(L*duration)/sum(duration) into :L_mean from simqueue;
quit;
	

data _null_;
	file print;
	title 'Simulation Results';
	put / @10 "Interarrival time distribution: &iatdistname (&iatparms)"
	    / @10 "Service time distribution: &servdistname (&servparms)"
		/ @10 "Number of servers: &nserv"
		/ @10 "Number of tasks: &ntask"
		// @10 "Interarrival time, mean: &iat_mean"
		/ @10 "Interarrival rate(lambda): %sysevalf(1/&iat_mean)"
		/ @10 "Service time, mean: &service_mean"
		/ @10 "Service rate (mu): %sysevalf(1/&service_mean)"
		/ @10 "Traffic intensity (ro): %sysevalf(&service_mean/(&iat_mean*&nserv))"
		/ @10 "Queue length, mean (Lq): &L_mean"
		/ @10 "Waiting time, mean (Wq): &wait_mean"
		/ @10 "Calculated queue length by Little's Law (lambda*Wq): %sysevalf(%sysevalf(1/&iat_mean)*&wait_mean)"
;
run;
title;

/*creating 1-row dataset with same stats*/
data outstat;
	nserv=&nserv;
	ntask=&ntask;
	lambda=%sysevalf(1/&iat_mean);
	service_mean=&service_mean;
	mu=%sysevalf(1/&service_mean);
	ro=%sysevalf(&service_mean/(&iat_mean*&nserv));
	IAT_est=&iat_mean;
	Ws_est=service_mean;
	Lq_est=&L_mean;
	Wq_est=&wait_mean;
	label
		nserv="Number of servers"
		ntask="Number of tasks"
		IAT_est="Interarrival time estimated, mean"
		lambda="Interarrival rate (lambda)"
		service_mean="Service time, mean"
		mu="Service rate (mu)"
		ro="Traffic intensity (ro)"
		Ws_est="Service time estimated, mean (Ws)"
		Lq_est="Queue length estimated, mean (Lq)"
		Wq_est="Waiting time estimated, mean (Wq)";
run;

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

options nomprint;

%mend qsim;


/*Example*/
%qsim(
  ntask=10000
 ,nserv=6
 ,iatdist=exponential(0.04)
 ,servdist=exponential(0.007)
 ,seed=-1
)