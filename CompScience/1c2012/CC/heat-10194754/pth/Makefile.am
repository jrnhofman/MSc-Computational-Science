noinst_PROGRAMS = heat 

heat_SOURCES = $(COMMON)/main.c compute.c

## Uncomment the following:
noinst_PROGRAMS += POSIX_barrier
POSIX_barrier_SOURCES = assignment_2.3/original.c

noinst_PROGRAMS += sem_barrier
sem_barrier_SOURCES = assignment_2.3/sema.c

noinst_PROGRAMS += mutex_barrier
mutex_barrier_SOURCES = assignment_2.3/mutex.c

noinst_PROGRAMS += sieve
sieve_SOURCES = assignment_2.4/sieve_fast.c

#
# Maintenance rules
#
CONFIGURE_DEPENDENCIES = $(srcdir)/../VERSION
ACLOCAL_AMFLAGS = --force -I .
include $(srcdir)/../src/libcommon.mk


