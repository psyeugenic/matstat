%%
%% File:    matstat.erl
%% Author:  Björn-Egil Dahlberg
%% Created: 2012-10-19
%%

-module(matstat).

-define(nolimit, inf).

-export([new/0, new/1, add/2]).

-export([msn/1,
         mean/1,
         divide/2,
         divide/3,
         mult/2,
         mult/3,
         cov/2,
         corr/2]).

-export([tmean/1, tmean/2,
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
         pearsonsr/1]).

-export([moment/2,
         kurtosis/1,
         skewness/1]).

-export([histogram/1, histogram/2,
         histogram_property/1]).

%%% stats implementation

-type limit() :: number() | 'inf'.
-type limit_min() :: limit().
-type limit_max() :: limit().

-record(moment, {
        m1 = 0,   % mean
        m2 = 0,
        m3 = 0,
        m4 = 0
    }).

-record(stats, {
        lmin   = ?nolimit :: limit_min(),
        lmax   = ?nolimit :: limit_max(),
        sum    = 0 :: number(),
        sumsq  = 0 :: number(),
        min    = 0 :: number(),
        max    = 0 :: number(),
        cbs    = [],
        n      = 0 :: non_neg_integer()
    }).

-type histogram() :: term().

-record(histogram, {
        n     = 0 :: non_neg_integer(),    % number of values in bins
        nbins = 0 :: non_neg_integer(),
        width = 1 :: non_neg_integer(),
        min   = none, % bin min
        max   = none, % bin max
        bins  = gb_trees:empty()
    }).


-type stats() :: term().

-spec new() -> stats().

new() -> new([]).

-spec new([Opt]) -> stats() when
      Opt :: {'min', limit_min()}
           | {'max', limit_max()}
           | {'histogram',number(),number(),number()}
           | 'gmean' | 'hmean'.

new(Opts) ->
    Cbs = [{moment, {fun update_moment/3, #moment{}}}],
    lists:foldl(fun(Opt,S) -> opts(Opt,S) end, #stats{cbs = Cbs}, Opts).

opts({min, V}, S) when is_number(V); V =:= ?nolimit -> S#stats{lmin = V};
opts({max, V}, S) when is_number(V); V =:= ?nolimit -> S#stats{lmax = V};
opts(gmean, #stats{cbs=Cbs}=S) -> S#stats{cbs=[{gmean, {fun update_gmean/3, 0}}|Cbs]};
opts(hmean, #stats{cbs=Cbs}=S) -> S#stats{cbs=[{hmean, {fun update_hmean/3, 0}}|Cbs]};
opts({histogram,Min,Max,N}, #stats{cbs=Cbs}=S) ->
    S#stats{cbs=[{histogram, {fun update_histogram/3, histogram_state(Min,Max,N)}}|Cbs]};
opts(_, S) -> S.

-spec add([number()] | number(), stats()) -> stats().

add(V, #stats{lmin = L, lmax = U} = S) when is_number(V) ->
    if ((L =:= ?nolimit orelse V >= L) andalso
        (U =:= ?nolimit orelse V =< U)) -> update(V, S);
       true -> S
    end;
add([V|Vs], S) -> add(Vs, add(V,S));
add(_, S) -> S.

update(V, #stats{n = 0, cbs = Cbs} = S) ->
    S#stats{
        min    = V,
        max    = V,
        n      = 1,
        sum    = V,
        sumsq  = V*V,
        cbs    = update_cbs(V, S, Cbs)
    };
update(V, #stats{ n = N, cbs = Cbs } = S) ->
    S#stats{
        min    = erlang:min(V, S#stats.min),
        max    = erlang:max(V, S#stats.max),
        n      = 1   + N,
        sum    = V   + S#stats.sum,
        sumsq  = V*V + S#stats.sumsq,
        cbs    = update_cbs(V, S, Cbs)
    }.

update_cbs(V, S, [{K,{Fun, Data}}|Cbs]) ->
    [{K, {Fun, Fun(V, S, Data)}}|update_cbs(V, S, Cbs)];
update_cbs(_, _, []) -> [].

update_moment(V, #stats{n = N1}, #moment{m1 = M1, m2 = M2, m3 = M3, m4 = M4} = Ms) ->
    N = N1 + 1,
    Delta = V - M1,
    DeltaN  = Delta / N,
    DeltaN2 = DeltaN*DeltaN,
    Term1 = Delta * DeltaN * N1,

    Ms#moment{
        m1 = M1 + DeltaN,
        m2 = M2 + Term1,
        m3 = M3 + Term1 * DeltaN * (N - 2) - 3 * DeltaN * M2,
        m4 = M4 + Term1 * DeltaN2 * (N*N - 3*N + 3) + 6 * DeltaN2 * M2 - 4 * DeltaN * M3
    }.

update_gmean(V,_,S) ->
    S + math:log(V).
update_hmean(V,_,S) ->
    S + 1 / V.


%%% vanilla implementation

-spec mean([number()] | stats()) -> float().

mean(Vs)  -> tmean(Vs).

-spec tmean([number()] | stats()) -> float().

tmean(#stats{n = N, sum = Sum}) -> Sum/N;
tmean(Vs) when is_list(Vs) -> tmean(Vs, {?nolimit, ?nolimit}).

-spec tmean([number()], {limit_min(),limit_max()}) -> float().

tmean(Vs, {L,U}) when is_list(Vs) ->
    tmean(add(Vs, new([{min,L},{max, U}]))).

%% Calculate nth root of (x1 * x2 * .. * xn)
-spec gmean([number()] | stats()) -> float().

gmean(#stats{n = N, cbs = Cbs}) when N > 0 ->
    {_, SL} = proplists:get_value(gmean, Cbs),
    math:exp(SL/N);
gmean(Vs) when is_list(Vs) -> gmean(Vs, {?nolimit,?nolimit}).

-spec gmean([number()], {limit_min(), limit_max()}) -> float().

gmean(Vs,{L,U}) when is_list(Vs) ->
    gmean(add(Vs, new([{min,L},{max,U},gmean]))).

%% Calculate n / (1/x1 + 1/x2 + .. + 1/xn)
%% x1 .. Xn > 0.0 (positive real numbers)

-spec hmean([number()]|stats()) -> float().

hmean(#stats{n=N,cbs=Cbs}) ->
    {_,S} = proplists:get_value(hmean,Cbs),
    N/S;
hmean(Vs) when is_list(Vs) -> hmean(Vs, {?nolimit, ?nolimit}).

-spec hmean([number()], {limit_min(), limit_max()}) -> float().

hmean(Vs, {L,U}) when is_list(Vs) ->
    hmean(add(Vs, new([{min,L},{max,U},hmean]))).


-spec tmin([number()] | stats()) -> number().

tmin(#stats{min = V}) -> V;
tmin(Vs) when is_list(Vs) -> tmin(Vs, ?nolimit).

-spec tmin([number()], limit_min()) -> number().

tmin(Vs, L) when is_list(Vs) ->
    tmin(add(Vs, new([{min, L}]))).

-spec tmax([number()] | stats()) -> number().

tmax(#stats{max = V }) -> V;
tmax(Vs) when is_list(Vs) -> tmax(Vs, ?nolimit).

-spec tmax([number()], limit_max()) -> number().

tmax(Vs, L) when is_list(Vs) ->
    tmax(add(Vs, new([{max, L}]))).

-spec tvar([number()] | stats()) -> float().

tvar(#stats{n = N, sum = Sum, sumsq = SumSqr}) ->
    (SumSqr - Sum*Sum/N)/(N - 1);
tvar(Vs) when is_list(Vs) -> tvar(Vs, {?nolimit, ?nolimit}).

-spec tvar([number()], {limit_min(), limit_max()}) -> float().

tvar(Vs, {L,U}) when is_list(Vs) ->
    tvar(add(Vs, new([{min,L},{max,U}]))).

-spec tstd([number()] | stats()) -> float().

tstd(#stats{} = S) -> math:sqrt(tvar(S));
tstd(Vs) when is_list(Vs) -> tstd(Vs, {?nolimit, ?nolimit}).

-spec tstd([number()], {limit_min(), limit_max()}) -> float().

tstd(Vs, {L,U}) when is_list(Vs) ->
    tstd(add(Vs, new([{min,L},{max,U}]))).

-spec tsem([number()] | stats()) -> float().

tsem(#stats{ n = N, sum = Sum, sumsq = SumSq }) ->
    math:sqrt(((SumSq - Sum*Sum/N)/(N - 1))/N);
tsem(Vs) when is_list(Vs) -> tsem(Vs, {?nolimit, ?nolimit}).

-spec tsem([number()], {limit_min(), limit_max()}) -> float().

tsem(Vs, {L,U}) when is_list(Vs) ->
    tsem(add(Vs, new([{min,L},{max,U}]))).

% from http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Higher-order_statistics

-spec kurtosis(stats() | [number()]) -> float().

kurtosis(#stats{ n = N } = S) when N > 0 ->
    M2 = moment(2, S),
    M4 = moment(4, S),
    G2 = M4 / (M2*M2) - 3, % true
    ((N + 1)*G2 + 6)*(N - 1)/((N-2)*(N-3)); % sampled
kurtosis(Vs) when is_list(Vs) -> kurtosis(Vs, {?nolimit, ?nolimit}).

-spec kurtosis([number()], {limit_min(), limit_max()}) -> float().

kurtosis(Vs, {L,U}) when is_list(Vs) ->
    kurtosis(add(Vs, new([{min,L},{max,U}]))).

-spec skewness(stats() | [number()]) -> float().

skewness(#stats{ n = N } = S) when N > 0 ->
    M2 = moment(2, S),
    M3 = moment(3, S),
    G1 = M3 / (M2 *math:sqrt(M2)), % true
    math:sqrt(N*(N-1))/(N - 2)*G1; % sampled
skewness(Vs) when is_list(Vs) -> skewness(Vs, {?nolimit, ?nolimit}).

-spec skewness([number()], {limit_min(), limit_max()}) -> float().

skewness(Vs, {L,U}) when is_list(Vs) ->
    skewness(add(Vs, new([{min,L},{max,U}]))).

-spec moment(Moment :: integer(), stats() | [number()]) -> float().

moment(I, #stats{ n = N, cbs = Cbs}) when is_integer(I), I > 0, I < 5->
    {_, M} = proplists:get_value(moment, Cbs),
    element(I + 1, M)/N;
moment(I, Vs) when is_list(Vs), is_integer(I) ->
    moment(I, add(Vs, new())).


%% histogram


-spec histogram([number()], integer()) -> [{number(), integer()}].

histogram(Vs, Nb) when is_integer(Nb) ->
    {Min, Max, _N}= min_max_n(Vs),
    histogram(add(Vs, new([{histogram,Min,Max,Nb}]))).

-spec histogram([number()]) -> [{number(), integer()}].

histogram([_|_]= Vs) ->
    {Min, Max, N}= min_max_n(Vs),
    Nb = if N > 4 -> trunc(math:sqrt(N)); true -> 2 end,
    histogram(add(Vs, new([{histogram,Min,Max,Nb}])));
histogram(#stats{cbs=Cbs}) ->
    {_,#histogram{nbins = N, min = Min,
               width = W, bins = Bins}} = proplists:get_value(histogram, Cbs),
    histogram_counts(0, N, Min, W, Bins).

histogram_counts(I, N, Min, W, Bins) when I < N ->
    [{I*W + Min, histogram_count(I, Bins)}|histogram_counts(I + 1, N, Min, W, Bins)];
histogram_counts(_, _, _, _, _) -> [].

histogram_count(Bin, Bins) ->
    case gb_trees:lookup(Bin, Bins) of
        none -> 0;
        {value, V} -> V
    end.

histogram_state(Min,Max,N) when is_number(Min), is_number(Max), is_integer(N), N > 1 ->
    Width = histogram_bin_width(Min,Max,N),
    #histogram{nbins = N, min = Min, max = Max, width = Width}.


update_histogram(V, _Stats, #histogram{n = N, bins = Bins0} = Hi) ->
    Ix = histogram_bin_idx(V, Hi),
    Bins1  = case gb_trees:lookup(Ix, Bins0) of
                 none       -> gb_trees:enter(Ix, 1, Bins0);
                 {value, C} -> gb_trees:enter(Ix, C + 1, Bins0)
             end,
    Hi#histogram{n = N + 1, bins = Bins1}.

histogram_bin_idx(V, #histogram{nbins = N, max = Max}) when V > Max -> N - 1;
histogram_bin_idx(V, #histogram{min = Min}) when V < Min -> 0;
histogram_bin_idx(V, #histogram{width = W, min = Min}) ->
    erlang:trunc((V - Min)/W).

-spec histogram_property(histogram()) ->
    [{atom(), number()}].

histogram_property(#histogram{n = N, width = W, min = Min, max = Max, nbins = Nbins}) ->
    [{n, N}, {width, W}, {min, Min}, {max, Max}, {nbins, Nbins}].

% try to keep width as integer if possible
histogram_bin_width(Min, Max, Nb) when is_integer(Min), is_integer(Max), is_integer(Nb) ->
    if (Max - Min + 1) rem Nb =:= 0 -> (Max - Min + 1) div Nb;
       true -> (Max - Min + 1)/Nb
    end;
histogram_bin_width(Min, Max, Nb) ->
    (Max - Min + 1)/Nb.


%% not in main framework

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
    RSQ   = (SSYY - SSE)/SSYY,
    SD    = if N > 2 -> math:sqrt(SSE/(N - 2));
              true -> 0.0
           end,
    {{Slope, Incpt}, {RSQ,SD}}.

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
            maybe_int_div(I0+I1,2)
    end.

maybe_int_div(T,N) when is_integer(T), is_integer(N) ->
    case T rem N of
        0 -> T div N;
        _ -> T / N
    end;
maybe_int_div(T,N) ->
    T / N.

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

min_max_n(Vs) ->
    V = first_value(Vs),
    min_max_n(Vs,{V,V,0}).

first_value(V) when is_number(V) -> V;
first_value([V|_]) -> first_value(V);
first_value(V) -> erlang:error({min_max_n,V}).

min_max_n([V|Vs], {Min, Max, N}) when is_number(V) ->
    Vmin = if V < Min -> V; true -> Min end,
    Vmax = if V > Max -> V; true -> Max end,
    min_max_n(Vs, {Vmin, Vmax, N + 1});
min_max_n([V|Vs], MMN) ->
    min_max_n(Vs,min_max_n(V, MMN));
min_max_n([], MMN) ->
    MMN.

sum_x_y_xy_x2_y2_n(Vs) -> 
    sum_x_y_xy_x2_y2_n(Vs, 0, 0, 0, 0, 0, 0).

sum_x_y_xy_x2_y2_n([{X,Y}|Vs], SumX, SumY, SumXY, SumX2, SumY2, N) when is_number(X), is_number(Y) ->
    sum_x_y_xy_x2_y2_n(Vs, SumX + X, SumY + Y, SumXY + X*Y, SumX2 + X*X, SumY2 + Y*Y, N + 1);
sum_x_y_xy_x2_y2_n([], SumX, SumY, SumXY, SumX2, SumY2, N) ->
    {SumX, SumY, SumXY, SumX2, SumY2, N}.


lists_tail([_|R], N) when N > 0 -> lists_tail(R, N - 1);
lists_tail(R, 0) -> R.
