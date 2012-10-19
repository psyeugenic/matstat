matstat
-------

Erlang statistics module. The current stdlib in OTP is lacking a statistics module.
I hope to add one sooner or later, in the meanwhile I will add stuff here.

In other words: the API will be subject to change until I feel it could be fixated.

Currently we have:

* `mean/1 -> Mean :: float()` - calculate the mean of list of numbers,
* `msn/1 -> {Mean :: float(), StdDev :: float()}` - calculate mean and *sampled* standard deviation,
* operations

Todo:

* histogram
* distributions
* anything useful for testing
