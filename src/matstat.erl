%%
%% Copyright (C) 2012 Hapida AB
%%
%% File:    matstat.erl
%% Author:  Björn-Egil Dahlberg
%% Created: 2012-10-19
%%

-module(matstat).
-export([
	msn/1,
	mean/1,
	divide/2,
	divide/3,
	mult/2,
	mult/3,
	cov/2,
	corr/2
    ]).


mean(Vs) -> mean(Vs, 0, 0).

mean([], N, Sum) -> Sum/N;
mean([V|Vs], N, Sum) -> mean(Vs, N + 1, Sum + V).


msn([])  -> {0, 0}; 
msn([V]) -> {V, 0.0};
msn(Vs)  -> msn(Vs, 0, 0, 0).

msn([V | Vs], Sum, SumSq, N) -> 
    msn(Vs, Sum + V, SumSq + V*V, N + 1);
msn([], Sum, SumSq, N) -> 
    Mean   = Sum / N,
    StdDev = math:sqrt((SumSq - (Sum*Sum/N))/(N - 1)),
    {Mean, StdDev}.


mult(A, B) -> mult(A, B, 0).
mult(A, B, C) -> divide(A, B, (-1.0)*C).

divide(A, B) -> divide(A, B, 0).
divide({Am, As}, {Bm, Bs}, C) ->
    F = Am/Bm,
    Fs2Fm2 = (As*As)/(Am*Am) + (Bs*Bs)/(Bm*Bm) - 2*(As*Bs)/(Am*Bm)*C,
    Fs = math:sqrt(Fs2Fm2*F*F),
    {F, Fs}.


corr(Xs, Ys) ->
   N  = length(Xs),
   {Xm, Xsd} = msn(Xs),
   {Ym, Ysd} = msn(Ys),
   C = cov(Xs, Xm, Ys, Ym, 0)/(N - 1),
   C/(Xsd*Ysd).

cov(Xs, Ys) ->
    N  = length(Xs),
    Xm = lists:sum(Xs)/N,
    Ym = lists:sum(Ys)/N,
    cov(Xs, Xm, Ys, Ym, 0)/(N - 1).

cov([], _, [], _, S) -> S;
cov([X|Xs], Xm, [Y|Ys], Ym, S) ->
    cov(Xs, Xm, Ys, Ym, S + (X - Xm)*(Y - Ym)).

