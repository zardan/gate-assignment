SETS
    I   'pucks, arc between arrival and departure node' /I1*I10/
    K   'gates' /K1*K3/
    S   'source node'   /S/
    T   'terminal node' /T/

ALIAS(I,J);

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

     t_dep(I)   'departing time of puck i'
     /I1 12.75
      I2 16.00
      I3 11.67
      I4 17.67
      I5 15.08
      I6 12.33
      I7 14.50
      I8  7.00
      I9  7.33
      I10 8.00/;

VARIABLES
    X(I,J,K)    'decision varaibles'
    obj;

BINARY VARIABLE X(I,J,K);

EQUATIONS
    puckConstr(I,J,K)       'only straight pucks'
    gateToPuck(I,J)         'have to assign gate to each puck'
    backwardFlow(I,J,K)
    objective;

    puckConstr(I,J,K)       $(ord(I) ne ord(J)) .. X(I,J,K) =e= 0;
    gateToPuck(I,J)         $(ord(I) eq ord(J)) .. sum(K, X(I,J,K)) =e= 1;
    backwardFlow(I,J,K)     $(t_arr(J) < t_dep(I)) .. X(I,J,K) =e= 0;
    objective               .. obj =e= 1;

MODEL basic /all/;

SOLVE basic USING MINLP MINIMIZING obj;

**DISPLAY totWork, totTime.l, x.l;


** noFlowToSource(I,S,K)   'no flox back to the source node'
** noFlowToSource(I,S,K)   .. X(I,S,K) =e= 0;

** backwardFlow(I,J,K)     $(t_dep(J) > t_arr(I)) .. X(J,I,K) =e= 0;
