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

-export([
	tmean/1, tmean/2,
	tmin/1, tmin/2,
	tmax/1, tmax/2,
	gmean/1,
	hmean/1
    ]).

-define(nolimit, inf).

-spec mean([number()]) -> float().
-spec tmean([number()]) -> float().
-spec tmean([number()], { number() | 'inf', number() | 'inf' }) -> float().

mean(Is) -> tmean(Is, {?nolimit, ?nolimit}).
tmean(Is) -> tmean(Is, {?nolimit, ?nolimit}).
tmean(Is, Limit) -> tmean(Is, Limit, 0, 0).
tmean([I|Is], {Ll, Ul} = Ls, S, N) when is_number(I),
					(Ll =:= ?nolimit orelse I >= Ll),
					(Ul =:= ?nolimit orelse I =< Ul) ->
    tmean(Is, Ls, S + I, N + 1);
tmean([_|Is], Ls, S, N) ->
    tmean(Is, Ls, S, N);
tmean([], _, 0, 0) -> 0.0;
tmean([], _, S, N) -> S/N.

-spec tmin([number()]) -> number().
-spec tmin([number()], 'inf' | number()) -> number().

tmin(Is) -> tmin(Is, ?nolimit).
tmin([I|Is], Limit) when I >= Limit; Limit =:= ?nolimit ->
    tmin(Is, Limit, I);
tmin([_|Is], Limit) -> 
    tmin(Is, Limit).

tmin([I|Is], Limit, Min) when I < Min andalso (I >= Limit orelse Limit =:= ?nolimit) ->
    tmin(Is, Limit, I);
tmin([_|Is], Limit, Min) ->
    tmin(Is, Limit, Min);
tmin([], _, Min) ->
    Min.

-spec tmax([number()]) -> number().
-spec tmax([number()], 'inf' | number()) -> number().

tmax(Is) -> tmax(Is, ?nolimit).
tmax([I|Is], Limit) when I =< Limit; Limit =:= ?nolimit ->
    tmax(Is, Limit, I);
tmax([_|Is], Limit) -> 
    tmax(Is, Limit).

tmax([I|Is], Limit, Max) when I > Max andalso (I =< Limit orelse Limit =:= ?nolimit) ->
    tmax(Is, Limit, I);
tmax([_|Is], Limit, Max) ->
    tmax(Is, Limit, Max);
tmax([], _, Max) ->
    Max.

-spec gmean([number()]) -> float().

%% Calculate nth root of (x1 * x2 * .. * xn)
gmean([I|Is]) when is_number(I) ->
    gmean(Is, I, 1).
gmean([I|Is], P, N) when is_number(I) ->
    gmean(Is, P*I, N + 1);
gmean([], P, N) ->
    math:pow(P, 1/N).


-spec hmean([number()]) -> float().

%% Calculate n / (1/x1 + 1/x2 + .. + 1/xn)
%% FIXME: zero numbers?
hmean(Is) -> hmean(Is, 0, 0).
hmean([I|Is], S, N) when is_number(I), I =/= 0 ->
    hmean(Is, S + (1/I), N + 1);
hmean([_|Is], S, N) ->
    hmean(Is, S, N);
hmean([], S, N) ->
    N / S.




%% old thinking
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

