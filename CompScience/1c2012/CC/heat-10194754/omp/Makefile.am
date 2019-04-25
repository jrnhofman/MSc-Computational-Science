noinst_PROGRAMS = heat

heat_SOURCES = $(COMMON)/main.c compute.c

noinst_PROGRAMS += matrix
matrix_SOURCES = matrix.c

#
# Maintenance rules
#
CONFIGURE_DEPENDENCIES = $(srcdir)/../VERSION
ACLOCAL_AMFLAGS = --force -I .
include $(srcdir)/../src/libcommon.mk


