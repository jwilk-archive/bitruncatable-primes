# Copyright © 2018-2019 Jakub Wilk <jwilk@jwilk.net>
# SPDX-License-Identifier: MIT

CXXFLAGS ?= -O2 -g
OPENMP_CXXFLAGS = -fopenmp
CXXFLAGS += -Wall $(OPENMP_CXXFLAGS)
LDLIBS += -lgmp -lgmpxx

.PHONY: all
all: bitruncatable-primes

.PHONY: clean
clean:
	rm -f bitruncatable-primes

# vim:ts=4 sts=4 sw=4 noet
