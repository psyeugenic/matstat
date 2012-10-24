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
	tvar/1, tvar/2,
	tstd/1, tstd/2,
	tsem/1, tsem/2,
	gmean/1,
	hmean/1,
	cmedian/1,
	linregress/1,
	itemfreq/1,
	pearsonsr/1
    ]).

-export([
	histogram/1, histogram/2,
	histogram_new/3,
	histogram_add/2,
	histogram_property/1,
	histogram_counts/1
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

-spec tvar([number()]) -> float().
-spec tvar([number()], {'inf' | number(), 'inf' | number()}) -> float().

tvar(Vs) -> tvar(Vs, {?nolimit, ?nolimit}).
tvar(Vs, Limit) ->
    {Sum,SumSqr,N} = sum_sumsqr_n(Vs, Limit, 0.0, 0.0, 0),
    (SumSqr - Sum*Sum/N)/(N - 1).

-spec tstd([number()]) -> float().
-spec tstd([number()], {'inf' | number(), 'inf' | number()}) -> float().

tstd(Vs) -> tstd(Vs, {?nolimit, ?nolimit}).
tstd(Vs, Limit) ->
    math:sqrt(tvar(Vs, Limit)).

-spec tsem([number()]) -> float().
-spec tsem([number()], {'inf' | number(), 'inf' | number()}) -> float().

tsem(Vs) -> tsem(Vs, {?nolimit, ?nolimit}).
tsem(Vs, Limit) ->
    {Sum,SumSqr,N} = sum_sumsqr_n(Vs, Limit, 0.0, 0.0, 0),
    math:sqrt(((SumSqr - Sum*Sum/N)/(N - 1))/N).
    


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
%% x1 .. Xn > 0.0 (positive real numbers)
hmean(Is) -> hmean(Is, 0, 0).
hmean([I|Is], S, N) when is_number(I), I > 0 ->
    hmean(Is, S + (1/I), N + 1);
hmean([_|Is], S, N) ->
    hmean(Is, S, N);
hmean([], S, N) ->
    N / S.


-spec linregress([{X :: number(), Y :: number()}]) ->
    {{Slope :: float(), Intercept :: float()}, {RSQ :: float(), StdDev :: float()}}.

%% Fix P-Value and what exact standard error really is?
linregress(Vs) ->
    {SumX, SumY, SumXY, SumX2, SumY2, N} = sum_x_y_xy_x2_y2_n(Vs),
    SSXY  = SumXY - SumX*SumY/N,
    SSXX  = SumX2 - SumX*SumX/N,
    Slope = SSXY/SSXX,
    Incpt = SumY/N - Slope*SumX/N,
    % error estimation
    SSYY  = SumY2 - SumY*SumY/N,
    SSE   = SSYY - Slope*SSXY,
    R2    = (SSYY - SSE)/SSYY,
    _SD    = if
	N > 2 -> math:sqrt(SSE/(N - 2));
	true -> 0.0
    end,
    {{Slope, Incpt}, R2}.

-spec pearsonsr([{X :: number(), Y :: number()}]) -> float().
    

pearsonsr(Vs) ->
    {SumX, SumY, SumXY, SumX2, SumY2, N} = sum_x_y_xy_x2_y2_n(Vs),
    DSxySxSy = (N *SumXY - SumX*SumY),
    Sqrx = math:sqrt(N*SumX2 - SumX*SumX),
    Sqry = math:sqrt(N*SumY2 - SumY*SumY),
    DSxySxSy/(Sqrx*Sqry).


-spec cmedian([number()]) -> number().

%% I think this should be implemented with histogram
%% could probably be done in k*O(n) instead of O(n * log n)

%% median
%% odd -> Xm
%% even -> (Xm_-1 + Xm_+1)/2
cmedian(Is) ->
    Ls = lists:sort(Is),
    N  = length(Ls),
    H  = N div 2,
    case N rem 2 of
	1 -> % odd
	    [I|_] = lists_tail(Ls, H),
	    I;
	0 -> % even
	    [I0,I1|_] = lists_tail(Ls, H - 1),
	    (I0 + I1) / 2
    end.

-spec itemfreq([term()]) -> [{term(), integer()}].

%% gb_trees uses coercion, e.g. 0 == 0.0
itemfreq(Vs) -> 
    lists:sort(fun
	    ({_,A}, {_,B}) when A > B -> true;
	    (_, _) -> false
	end, gb_trees:to_list(itemfreq(Vs, gb_trees:empty()))).

itemfreq([V|Vs], T) ->
    case gb_trees:lookup(V, T) of
	none -> itemfreq(Vs, gb_trees:enter(V, 1, T));
	{value, N} -> itemfreq(Vs, gb_trees:enter(V, N + 1, T))
    end;
itemfreq([], T) -> T.



-record(hist, {
	n     = 0,    % number of values in bins
	nbins = 0,
	width = 1,
	min   = none, % bin min
	max   = none, % bin max
	bins  = gb_trees:empty()
    }).

-type histogram() :: term().

-spec histogram([number()]) -> [{number(), integer()}].

histogram([_|_]= Vs) ->
    {CVs, Min, Max, N}= min_max_n(Vs),
    Nb = if
	N > 4 -> trunc(math:sqrt(N));
	true -> 2
    end,
    histogram_counts(histogram_add(CVs, histogram_new(Min, Max, Nb))).

-spec histogram([number()], integer()) -> [{number(), integer()}].

histogram(Vs, Nb) when is_integer(Nb) ->
    {CVs, Min, Max, _N}= min_max_n(Vs),
    histogram_counts(histogram_add(CVs, histogram_new(Min, Max, Nb))).

-spec histogram_new(number(), number(), integer()) -> histogram().

histogram_new(Min, Max, Nb) when is_number(Min), is_number(Max), is_integer(Nb), Nb > 1 ->
    Width = histogram_bin_width(Min, Max, Nb),
    #hist{ nbins = Nb, min = Min, max = Max, width = Width }.

-spec histogram_add(number() | [number()], histogram()) -> histogram().

histogram_add(V, #hist{ n = N, bins = Bins0 } = Hi) when is_number(V) ->
    BinIdx = histogram_bin_idx(V, Hi),
    Bins1  = case gb_trees:lookup(BinIdx, Bins0) of
	none       -> gb_trees:enter(BinIdx, 1, Bins0);
	{value, C} -> gb_trees:enter(BinIdx, C + 1, Bins0)
    end,
    Hi#hist{ n = N + 1, bins = Bins1 };

histogram_add([], Hs) -> Hs;
histogram_add([V|Vs], Hs) ->
    histogram_add(Vs, histogram_add(V, Hs)).

-spec histogram_property(histogram()) ->
    [{atom(), number()}].

histogram_property(#hist{ n = N, width = W, min = Min, max = Max, nbins = Nbins}) ->
    [ {n, N}, {width, W}, {min, Min}, {max, Max}, {nbins, Nbins} ].

-spec histogram_counts(histogram()) -> [{number(), integer()}].

histogram_counts(#hist{ nbins = N, min = Min, width = W, bins = Bins }) ->
    histogram_counts(0, N, Min, W, Bins).

histogram_counts(I, N, Min, W, Bins) when I < N ->
    [{I*W + Min, histogram_count(I, Bins)}|histogram_counts(I + 1, N, Min, W, Bins)];
histogram_counts(_, _, _, _, _) -> [].

histogram_count(Bin, Bins) ->
    case gb_trees:lookup(Bin, Bins) of
	none -> 0;
	{value, V} -> V
    end.

histogram_bin_idx(V, #hist{ nbins = N, max = Max }) when V > Max -> N - 1;
histogram_bin_idx(V, #hist{ min = Min }) when V < Min -> 0;
histogram_bin_idx(V, #hist{ width = W, min = Min }) ->
    erlang:trunc((V - Min)/W).

% try to keep width as integer if possible
histogram_bin_width(Min, Max, Nb) when is_integer(Min), is_integer(Max), is_integer(Nb) ->
    if
	(Max - Min + 1) rem Nb =:= 0 -> (Max - Min + 1) div Nb;
	true -> (Max - Min + 1)/Nb
    end;
histogram_bin_width(Min, Max, Nb) ->
    (Max - Min + 1)/Nb.


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


%% aux

min_max_n([V|Vs]) when is_number(V) -> min_max_n(Vs,[V],V,V,1);
min_max_n([_|Vs]) -> min_max_n(Vs).

min_max_n([V|Vs], Os, Min, Max, N) when is_number(V) ->
    Vmin = if V < Min -> V; true -> Min end,
    Vmax = if V > Max -> V; true -> Max end,
    min_max_n(Vs, [V|Os], Vmin, Vmax, N + 1);
min_max_n([_|Vs], Os, Min, Max, N) ->
    min_max_n(Vs, Os, Min, Max, N);
min_max_n([], Os, Min, Max, N) ->
    {Os, Min, Max, N}.

sum_x_y_xy_x2_y2_n(Vs) -> 
    sum_x_y_xy_x2_y2_n(Vs, 0, 0, 0, 0, 0, 0).

sum_x_y_xy_x2_y2_n([{X,Y}|Vs], SumX, SumY, SumXY, SumX2, SumY2, N) when is_number(X), is_number(Y) ->
    sum_x_y_xy_x2_y2_n(Vs, SumX + X, SumY + Y, SumXY + X*Y, SumX2 + X*X, SumY2 + Y*Y, N + 1);
sum_x_y_xy_x2_y2_n([], SumX, SumY, SumXY, SumX2, SumY2, N) ->
    {SumX, SumY, SumXY, SumX2, SumY2, N}.


sum_sumsqr_n([V|Vs], {Ll, Ul} = Ls, Sum, SumSqr, N) when is_number(V),
					(Ll =:= ?nolimit orelse V >= Ll),
					(Ul =:= ?nolimit orelse V =< Ul) ->
    sum_sumsqr_n(Vs, Ls, Sum + V , SumSqr + V*V, N + 1);
sum_sumsqr_n([_|Vs], Ls, Sum, SumSqr, N) ->
    sum_sumsqr_n(Vs, Ls, Sum, SumSqr, N);
sum_sumsqr_n([], _, Sum, SumSqr, N) ->
    {Sum, SumSqr, N}.



lists_tail([_|R], N) when N > 0 -> lists_tail(R, N - 1);
lists_tail(R, 0) -> R.
