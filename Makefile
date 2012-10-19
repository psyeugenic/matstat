ERL ?= erl
APP := matstat

.PHONY: deps test

all: deps
	@./rebar compile

deps:
	@./rebar get-deps

clean:
	@./rebar clean

distclean: clean
	@./rebar delete-deps

test:
	@./rebar ct

docs:
	@erl -noshell -run edoc_run application '$(APP)' '"."' '[]'
