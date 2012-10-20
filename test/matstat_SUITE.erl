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
-export([
	mean_stddev/1,
	mean/1,
	tmin/1,tmax/1,
	tstd/1, tvar/1,
	hmean/1,
	gmean/1,
	cmedian/1
    ]).

init_per_suite(Config) when is_list(Config) ->
    Config.

end_per_suite(Config) when is_list(Config) ->
    ok.

all() ->
    [
	mean, mean_stddev,
	hmean, gmean,
	tmin,tmax,
	tvar, tstd,
	cmedian
].

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
    M   = matstat:gmean(Set),
    ok  = ?equal(1.442249, M).

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
    ok.


cmedian(_Config) ->
    Set1 = [1,2,3,4,5],
    ok   = ?equal(3, matstat:cmedian(Set1)),
    Set2 = [2, -1.0, -1.1, 0, 1, 2],
    ok   = ?equal(0.5, matstat:cmedian(Set2)),
    Set3 = [10,-9,8,-3,1,3,4,5,3],
    ok   = ?equal(3, matstat:cmedian(Set3)),
    ok.

