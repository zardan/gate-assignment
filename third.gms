SETS
    I   arrival nodes /I1*I10/
    J   'departure nodes' /J1*J10/
    K   'gates' /K1*K6/
    S   'source node'   /S/
    T   'terminal node' /T/

    ALIAS(I,I2)
    ALIAS(J,J2)
    ALIAS(K,K2);

SCALAR
    buffer      'buffer time for gates to take new pucks'
        /0.25/
    gate_cost   'cost for using one gate for one day'
        /100/
    walk_cost   'cost for walking at a higher speed'
        /50/
    alpha       'weight on multiobjective function'
        /0/;

PARAMETERS
    t_arr(I)    'arrival time of puck I'
        /I1 11.88
         I2 15.05
         I3 10.96
         I4 16.85
         I5 14.23
         I6 11.55
         I7 13.43
         I8  6.50
         I9  6.83
         I10 7.50/

     t_dep(J)   'departing time of puck J'
        /J1 12.75
         J2 16.00
         J3 11.67
         J4 17.67
         J5 15.08
         J6 12.33
         J7 14.50
         J8  7.00
         J9  7.33
         J10 8.00/

      walk_dist(K,K2)  'cost of walking between gates'
      /K1 .K2 10,   K1 .K3 20,  K1 .K4 22,  K1 .K5 40,  K1 .K6 41
                    K2 .K3 22,  K2 .K4 20,  K2 .K5 41,  K2 .K6 40
                                K3 .K4 10,  K3 .K5 20,  K3 .K6 22
                                            K4 .K5 22,  K4 .K6 20
                                                        K5 .K6 10/

      no_conn_pass(I,I2) 'number of connecting passengers between flight I and flight J'
      /I3 .I5 10, I3 .I2 20
       I7 .I2 20, I7 .I4 10, I7 .I5 10/

      compatibility(I,K) 'imposing physical restrictions on flights and gates'
      /I1 .K1 1, I1 .K2 1, I1 .K3 1
       I2 .K1 1, I2 .K2 1, I2 .K3 1
       I3 .K4 1, I3 .K5 1, I3 .K6 1
       I4 .K1 1, I4 .K2 1, I4 .K3 1
       I5 .K4 1, I5 .K5 1, I5 .K6 1
       I6 .K4 1, I6 .K5 1, I6 .K6 1
       I7 .K1 1, I7 .K2 1, I7 .K3 1
       I8 .K4 1, I8 .K5 1, I8 .K6 1
       I9 .K4 1, I9 .K5 1, I9 .K6 1
       I10 .K1 1, I10 .K2 1, I10 .K3 1/

      t_diff(I,J)       'time between flights';
      t_diff(I,J) = t_dep(J)-t_arr(I);

      walk_dist(K,K2)$(ord(K) gt ord (K2)) = walk_dist(K2,K);


VARIABLES
    Xs(S,I,K)       'dec. var. from source node to arrival node'
    Xf(I,J,K)       'dec. var. from arrival node to departure node'
    Xb(J,I,K)       'dec. var. from departure node to arrival node'
    Xt(J,T,K)       'dec. var. from departure node to termnal node'
    Xst(S,T,K)      'dec. var. from source node to terminal node'
    Ind(J,I)        'indicator function. true (1) if possible with flow from J to I'

    gateUtilization     'cost assosiated with gate utilization'
    unconvenienceCost   'cost associates with unconvenience for connecting passengers'
    objective           'total objective';

BINARY VARIABLE Xs(S,I,K), Xf(I,J,K), Xb(J,I,K), Xt(J,T,K), Xst(S,T,K);

EQUATIONS
    sourceOutflow(S,T,K)        'all gates either go directly to terminal or to arrival node'
    sourceToArr(S,I)            'at most one gate to each arrival node'
    cofArr(I,S,K)               'conservation of flow (cof) at arrival node'
    cofDep(J,T,K)               'conservation of flow (cof) at departure node'
    bfIndicatorTrue(I,J)        'indicator function for backward flow. 1 if possible'
    bfIndicatorFalse(I,J)       'indicator function for backward flow. 0 if not possible'
    forwardFlowTrue(I,J)        'one gate for each puck'
    forwardFlowFalse(I,J)       'pucks are always straight (I=J)'
    backwardFlowTrue(I,J,K)     'backward flow only possible when indicator function nonzero'
    backwardFlowLim(I,J)        'at most one gate per backward flow'
    terminalInflow(S,T,K)       'all gates must come back to terminal'

    objectiveEq                 'total objective'
    gateUtilizationEq           'cost assosiated with gate utilization'
    unconvenienceCostEq         'cost associated with unconvenience for connecting passengers';

    sourceOutflow(S,T,K)        .. sum(I, compatibility(I,K)*Xs(S,I,K)) + Xst(S,T,K) =e= 1;
    sourceToArr(S,I)            .. sum(K, compatibility(I,K)*Xs(S,I,K)) =l= 1;
    cofArr(I,S,K)               .. compatibility(I, K)*Xs(S,I,K) + sum(J, Ind(J,I)*compatibility(I,K)*Xb(J,I,K)) =e= sum(J, Xf(I,J,K));
    cofDep(J,T,K)               .. Xt(J,T,K) + sum(I, Ind(J,I)*compatibility(I,K)*Xb(J,I,K)) =e= sum(I, Xf(I,J,K));
    bfIndicatorTrue(I,J)        $(t_arr(I)-t_dep(J)>buffer) .. Ind(J,I) =e= 1;
    bfIndicatorFalse(I,J)       $(t_arr(I)-t_dep(J)<buffer) .. Ind(J,I) =e= 0;
    forwardFlowTrue(I,J)        $(ord(I) eq ord(J)).. sum(K, Xf(I,J,K)) =e= 1;
    forwardFlowFalse(I,J)       $(ord(I) ne ord(J)).. sum(K, Xf(I,J,K)) =e= 0;
    backwardFlowTrue(I,J,K)     .. compatibility(I,K)*Xb(J,I,K) =l= Ind(J,I);
    backwardFlowLim(I,J)        .. sum(K, compatibility(I,K)*Xb(J,I,K)) =l= 1;
    terminalInflow(S,T,K)       .. sum(J, Xt(J,T,K)) + Xst(S,T,K) =e= 1;

    objectiveEq                 .. objective =e=  (1-alpha)*gateUtilization + alpha*unconvenienceCost;
    gateUtilizationEq           .. gateUtilization =e= sum((S,I,K), gate_cost*Xs(S,I,K));
    unconvenienceCostEq         .. unconvenienceCost =e= sum((I,I2,J,J2,K,K2), walk_cost*(walk_dist(K,K2)/t_diff(I,J2))*no_conn_pass(I,I2)*Xf(I,J,K)*Xf(I2,J2,K2));

MODEL basic /all/;

SOLVE basic USING MINLP MINIMIZING objective;