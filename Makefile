#
# Created by gmakemake (Ubuntu Mar 10 2012) on Thu Jul 26 14:42:21 2012
#
# Edited by arlange	
#
# Definitions
#

.SUFFIXES:
.SUFFIXES:	.a .o .c .C .cpp .s .S
.c.o:
		$(COMPILE.c) $<
.C.o:
		$(COMPILE.cc) $<
.cpp.o:
		$(COMPILE.cc) $<
.S.s:
		$(CPP) -o $*.s $<
.s.o:
		$(COMPILE.s) -o $@ $<
.c.a:
		$(COMPILE.c) -o $% $<
		$(AR) $(ARFLAGS) $@ $%
		$(RM) $%
.C.a:
		$(COMPILE.cc) -o $% $<
		$(AR) $(ARFLAGS) $@ $%
		$(RM) $%
.cpp.a:
		$(COMPILE.cc) -o $% $<
		$(AR) $(ARFLAGS) $@ $%
		$(RM) $%

AS =		as
CC =		gcc
CXX =		g++

RM = rm -f
AR = ar
LINK.c = $(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS)
LINK.cc = $(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) 
COMPILE.s = $(AS) $(ASFLAGS)
COMPILE.c = $(CC) $(CFLAGS) $(CPPFLAGS) -c
COMPILE.cc = $(CXX) $(CXXFLAGS) $(CPPFLAGS) -c
CPP = $(CPP) $(CPPFLAGS)
########## Default flags (redefine these with a header.mak file if desired)
CXXFLAGS =	-ggdb #-pg #flag for gprof
CFLAGS =	-ggdb
CLIBFLAGS =	-lm
CCLIBFLAGS =	
########## End of default flags


CPP_FILES =	g_core.cpp set_operations2.cpp gluier.cpp g_is.cp g_ir.cpp
C_FILES =	
PS_FILES =	
S_FILES =	
H_FILES =	g.h set_operations2.h
SOURCEFILES =	$(H_FILES) $(CPP_FILES) $(C_FILES) $(S_FILES)
.PRECIOUS:	$(SOURCEFILES)
OBJFILES =	g_core.o g_is.o g_ir.o set_operations2.o 
nautpath = 	/home/alex/software/nauty24r2/
NAUTOBJ = 	/home/alex/software/nauty24r2/gtools.o 
NAUTOBJ2 = 	$(nautpath)/gtnauty.o $(nautpath)/nauty.o $(nautpath)/nautil.o $(nautpath)/nautinv.o  $(nautpath)/naugraph.o $(nautpath)/nausparse.o

#
# Main targets
#

all:	gluier

gluier:	gluier.o $(OBJFILES)
	$(CXX) $(CXXFLAGS) -o gluier gluier.o $(OBJFILES) $(NAUTOBJ) $(NAUTOBJ2) $(CCLIBFLAGS)



#
# Dependencies
#

gluier.o:	g.h set_operations2.h

#
# Housekeeping
#

Archive:	archive.tgz

archive.tgz:	$(SOURCEFILES) Makefile
	tar cf - $(SOURCEFILES) Makefile | gzip > archive.tgz

clean:
	-/bin/rm -f $(OBJFILES) gluier.o *.o

realclean:        clean
	-/bin/rm -f gluier
