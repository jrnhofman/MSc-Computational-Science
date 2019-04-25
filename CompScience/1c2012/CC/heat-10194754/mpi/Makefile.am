noinst_PROGRAMS = heat

heat_SOURCES = $(COMMON)/main.c compute.c

## Uncomment the following:
noinst_PROGRAMS += benchmark
benchmark_SOURCES = benchmark.c

noinst_PROGRAMS += broadcast_test
broadcast_test_SOURCES = broadcast_test.c

#
# Maintenance rules
#
CONFIGURE_DEPENDENCIES = $(srcdir)/../VERSION
ACLOCAL_AMFLAGS = --force -I .
include $(srcdir)/../src/libcommon.mk


