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
	mean/1,
	mean_stddev/1,
	gmean/1
    ]).

init_per_suite(Config) when is_list(Config) ->
    Config.

end_per_suite(Config) when is_list(Config) ->
    ok.


all() ->
    [mean, mean_stddev,
    gmean].

-define(err, (0.00005)).
-define(equal(A,B), if erlang:abs(A-B) < ?err -> ok; true -> {A,B} end).

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
    ok    = ?equal(19.1311, S).

gmean(_Config) ->
    Set = [1.8, 1.16666, 1.428571],
    M   = matstat:gmean(Set),
    ok  = ?equal(1.442249, M).

