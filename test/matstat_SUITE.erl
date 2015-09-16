%%
%% Copyright (C) 2012 Hapida AB
%%
%% File:    matstat_SUITE.erl
%% Author:  Björn-Egil Dahlberg
%% Created: 2012-10-19
%%

-module(matstat_SUITE).
-include("test_server.hrl").

%% common test
-export([all/0]).
-export([init_per_suite/1, end_per_suite/1]).

%% Test cases
-export([mean_stddev/1,
        mean/1,
        tmin/1,tmax/1,
        tstd/1, tvar/1,
        tsem/1,
        hmean/1,
        gmean/1,
        cmedian/1,
        linregress/1,
        itemfreq/1,
        pearsonsr/1,
        moment/1,
        skewness/1,
        kurtosis/1,
        histogram/1]).

init_per_suite(Config) when is_list(Config) ->
    Config.

end_per_suite(Config) when is_list(Config) ->
    ok.

all() -> [mean, mean_stddev,
          hmean, gmean,
          tmin,tmax,
          tvar, tstd,
          tsem,
          cmedian,
          linregress,
          itemfreq,
          pearsonsr,
          moment,
          skewness,
          kurtosis,
          histogram].

-define(err, (0.005)).
-define(equal(A,B), equal(A,B)).

equal(A,B) when abs(A-B) < ?err -> ok;
equal(A,B) -> {A,B}.

%%----------------------------------------------------------------------
%% Tests
%%----------------------------------------------------------------------

mean(suite) -> [];
mean(_Config) ->
    Set  = [70, 82, 76, 79, 83, 85, 72, 77],
    M    = matstat:mean(Set),
    ok   = ?equal(78, M).

mean_stddev(suite) -> [];
mean_stddev(_Config) ->
    Set   = [96, 104, 126, 134, 140],
    {M,S} = matstat:msn(Set),
    ok    = ?equal(120, M),
    ok    = ?equal(19.1311, S),
    ok.

gmean(_Config) ->
    Set = [1.8, 1.16666, 1.428571],
    M1  = matstat:gmean(Set),
    ok  = ?equal(1.442249, M1),
    S   = lists:foldl(fun(V,Si) -> matstat:add(V,Si) end, matstat:new([gmean]), Set),
    M2  = matstat:gmean(S),
    ok  = ?equal(1.442249, M2),
    ok.

hmean(_Config) ->
    Set = [1,2,4],
    M   = matstat:hmean(Set),
    ok  = ?equal(12/7, M).

tmin(_Config) ->
    Set1 = [1,2,3,4,5],
    ok   = ?equal(1, matstat:tmin(Set1)),
    ok   = ?equal(3, matstat:tmin(Set1, 3)),
    Set2 = [2, -1.0, -1.1, 0, 1, 2],
    ok   = ?equal(-1.0, matstat:tmin(Set2, -1)),
    Set3 = [10,-9,8,-3,1,3,4,5,3],
    ok   = ?equal(-9, matstat:tmin(Set3)),
    ok   = ?equal(-3, matstat:tmin(Set3, -5)),
    ok.

tmax(_Config) ->
    Set1 = [1,2,3,4,5],
    ok   = ?equal(5, matstat:tmax(Set1)),
    ok   = ?equal(3, matstat:tmax(Set1,3)),
    Set2 = [2, -1.0, -1.1, 0, 1, 2],
    ok   = ?equal(2, matstat:tmax(Set2)),
    Set3 = [10,-9,8,-3,1,3,4,5,3],
    ok   = ?equal(10, matstat:tmax(Set3)),
    ok   = ?equal(-3, matstat:tmax(Set3, -2)),
    ok.

tvar(_Config) ->
    Set1 = lists:flatten([[1,2,3,4,5,6] || _ <- [1,2,3,4,5]]),
    ok   = ?equal(3.0172, matstat:tvar(Set1)),
    Set2 = [0,-1,-3,7,8,9] ++ Set1,
    ok   = ?equal(3.0172, matstat:tvar(Set2, {1,6})),
    ok.


tstd(_Config) ->
    Set1 = [73.0, 93.3, 183.7, 86.6, 77.3],
    ok   = ?equal(45.92, matstat:tstd(Set1)),
    ok   = ?equal(45.92, matstat:tstd(Set1,{inf,inf})),
    Set2 = [33,36,73.0, 93.3, 183.7, 86.6, 77.3, 203, 204],
    ok   = ?equal(45.92, matstat:tstd(Set2,{40,200})),
    Set3 = [9,2,1],
    ok   = ?equal(4.36, matstat:tstd(Set3)),
    Set4 = [6,6,8,8,3,8],
    ok   = ?equal(1.97, matstat:tstd(Set4)),
    ok.

tsem(_Config) ->
    Set1 = [9,2,1],
    ok   = ?equal(2.52, matstat:tsem(Set1)),
    Set2 = [6,6,8,8,3,8],
    ok   = ?equal(0.81, matstat:tsem(Set2)),
    Set3 = [1,2,1,10,11,-1,-1.0] ++ Set2,
    ok   = ?equal(0.81, matstat:tsem(Set3, {3,9})),
    ok.

cmedian(_Config) ->
    Set1 = [1,2,3,4,5],
    ok   = ?equal(3, matstat:cmedian(Set1)),
    Set2 = [2, -1.0, -1.1, 0, 1, 2],
    ok   = ?equal(0.5, matstat:cmedian(Set2)),
    Set3 = [10,-9,8,-3,1,3,4,5,3],
    ok   = ?equal(3, matstat:cmedian(Set3)),
    ok.

linregress(_Config) ->
    Set1 = [{95,214},{82,152},{90,156},{81,129},{99,254},{100,266},
            {93,210},{95,204},{93,213},{87,150}],
    {{Slope, Intercept}, {R2,SD}} = matstat:linregress(Set1),
    io:format("StdDev ~p~n", [SD]),
    ok   = ?equal(6.7175, Slope),
    ok   = ?equal(-419.85, Intercept),
    ok   = ?equal(0.89, R2),
    ok.

pearsonsr(_Config) ->
    Set = [{95,214},{82,152},{90,156},{81,129},{99,254},{100,266},
	{93,210},{95,204},{93,213},{87,150}],
    R  = matstat:pearsonsr(Set),
    ok = ?equal(0.944093033, R),
    ok.

itemfreq(_Config) ->
    Set1 = [1,1,1,1,1,a,a,a,a,d,d,d,e,e,"hi"],
    [{1,5},{a,4},{d,3},{e,2},{"hi",1}] = matstat:itemfreq(Set1),
    Set2 = [1,a,1,e,1,a,a,1,a,d,"hi",d,1,e,d,1,1,1],
    [{1,8},{a,4},{d,3},{e,2},{"hi",1}] = matstat:itemfreq(Set2),
    [] = matstat:itemfreq([]),
    ok.


skewness(_Config) ->
    % sampled
    Set1 = [m(61,5),m(64,18),m(67,42),m(70,27),m(73,8)],
    ok   = ?equal(-0.1098, matstat:skewness(Set1)),
    ok.

moment(_Config) ->
    Set1 = [1,3,6,10],
    Set2 = [m(61,5),m(64,18),m(67,42),m(70,27),m(73,8)],
    ok   = ?equal(8.5275, matstat:moment(2, Set2)),
    ok   = ?equal(-2.6933, matstat:moment(3, Set2)),
    ok   = ?equal(199.3760, matstat:moment(4, Set2)),

    % was these wrong to begin with?
    %ok   = ?equal(5, matstat:moment(1, Set1)),
    %ok   = ?equal(36.5, matstat:moment(2, Set1)),
    %ok   = ?equal(311, matstat:moment(3, Set1)),
    ok.

%moment_mean(_Config) ->
%    Set1 = [1,3,6,10],
%    ok   = ?equal(11.5, matstat:moment(2, mean, Set1)),
%    ok.

kurtosis(_Config) ->
    Set1 = [6,7,8,9,10],
    Set2 = [m(61,5),m(64,18),m(67,42),m(70,27),m(73,8)],
    % total, i.e. not sampled values
    % Set3 = [2, 2, 4, 6, 8, 10, 10],
    % ok   = ?equal(-(85/54), matstat:kurtosis(Set3)),
    % Set4 = [0, 7, 7, 6, 6, 6, 5, 5, 4, 1],
    % ok   = ?equal(-(74146/271441), matstat:kurtosis(Set4)),
    % sampled
    ok   = ?equal(-1.2, matstat:kurtosis(Set1)),
    ok   = ?equal(-0.2091, matstat:kurtosis(Set2)),
    ok.


histogram(_Config) ->
    Set1 = lists:seq(1,10),
    Set2 = [m(61,5),m(64,18),m(67,42),m(70,27),m(73,8)],
    _ = matstat:histogram(Set1),
    _ = matstat:histogram(Set2),
    Set3 = lists:seq(11,20),
    Set4 = lists:seq(21,30) ++ Set1 ++ Set3,
    S1   = matstat:new([{histogram,11,20,3}]),
    S2   = matstat:add(Set3,S1),
    _    = matstat:histogram(S2),
    S3   = matstat:add(Set4,S2),
    _    = matstat:histogram(S3),
    ok.


m(V,N) -> lists:duplicate(N,V).
