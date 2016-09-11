SETS
    I   'arrival nodes' /I1*I10/
    J   'departure nodes' /J1*J10/
    K   'gates' /K1*K5/
    S   'source node'   /S/
    T   'terminal node' /T/;

SCALAR
    buffer /0.25/
    gateCost /10/;

PARAMETERS
    t_arr(I)    'arrival time of puck i'
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

     t_dep(J)   'departing time of puck i'
     /J1 12.75
      J2 16.00
      J3 11.67
      J4 17.67
      J5 15.08
      J6 12.33
      J7 14.50
      J8  7.00
      J9  7.33
      J10 8.00/;



VARIABLES
    Xs(S,I,K)   'dec. var. from source node to arrival node'
    Xf(I,J,K)   'dec. var. from arrival node to departure node'
    Xb(J,I,K)   'dec. var. from departure node to arrival node'
    Xt(J,T,K)   'dec. var. from departure node to termnal node'
    Xst(S,T,K)  'dec. var. from source node to terminal node'
    Indic(J,I)    'indicatror!!'
    obj;

BINARY VARIABLE Xs(S,I,K), Xf(I,J,K), Xb(J,I,K), Xt(J,T,K), Xst(S,T,K);


EQUATIONS
    objective                   ''

    flowFromSource(S,T,K)       ''
    fromSourceToArr(S,I)        ''
    indicator1(I,J)             ''
    indicator2(I,J)             ''
    flowArrivalNode(I,S,K)      ''
    backwardFlowTrue(I,J,K)     ''
    flowDepartureNode(J,T,K)    ''
    oneGatePerXb(I,J)           ''
    forwardFlowTrue(I,J)        ''
    forwardFlowFalse(I,J)       ''
    flowToTerminal(S,T,K)       '';


    objective                   .. obj =e= sum((I,K,S), Xs(S,I,K))*gateCost;

    flowFromSource(S,T,K)       .. sum(I,Xs(S,I,K)) + Xst(S,T,K) =e= 1;
    fromSourceToArr(S,I)        .. sum(K, Xs(S,I,K)) =l= 1;


    indicator1(I,J)             $(t_arr(I)-t_dep(J)>buffer) .. Indic(J,I) =e= 1;
    indicator2(I,J)             $(t_arr(I)-t_dep(J)<buffer) .. Indic(J,I) =e= 0;
    flowArrivalNode(I,S,K)      .. Xs(S,I,K) + sum(J, Indic(J,I)*Xb(J,I,K)) =e= sum(J, Xf(I,J,K));
    backwardFlowTrue(I,J,K)     .. Xb(J,I,K) =l= Indic(J,I);
    flowDepartureNode(J,T,K)    .. Xt(J,T,K) + sum(I, Indic(J,I)*Xb(J,I,K)) =e= sum(I, Xf(I,J,K));
    oneGatePerXb(I,J)           .. sum(K, Xb(J,I,K)) =l= 1;
    forwardFlowTrue(I,J)        $(ord(I) eq ord(J)).. sum(K, Xf(I,J,K)) =e= 1;
    forwardFlowFalse(I,J)       $(ord(I) ne ord(J)).. sum(K, Xf(I,J,K)) =e= 0;
    flowToTerminal(S,T,K)       .. sum(J, Xt(J,T,K)) + Xst(S,T,K) =e= 1;


MODEL basic /all/;

SOLVE basic USING MINLP MINIMIZING obj;

**
