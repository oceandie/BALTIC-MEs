# Automatic Makefile

FCM_BLD_CFG = /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/mk/bld_tools.cfg
export FCM_VERBOSE ?= 1

export FCM_ROOTDIR := /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/BLD
export FCM_CACHEDIR := $(FCM_ROOTDIR)/.cache/.bld
export FCM_CFGDIR := $(FCM_ROOTDIR)/cfg
export FCM_SRCDIR := $(FCM_ROOTDIR)/src
export FCM_BINDIR := $(FCM_ROOTDIR)/bin
export FCM_ETCDIR := $(FCM_ROOTDIR)/etc
export FCM_LIBDIR := $(FCM_ROOTDIR)/lib
export FCM_TMPDIR := $(FCM_ROOTDIR)/tmp
export FCM_DONEDIR := $(FCM_ROOTDIR)/done
export FCM_FLAGSDIR := $(FCM_ROOTDIR)/flags
export FCM_INCDIR := $(FCM_ROOTDIR)/inc
export FCM_PPSRCDIR := $(FCM_ROOTDIR)/ppsrc
export FCM_OBJDIR := $(FCM_ROOTDIR)/obj
export FCM_ROOTPATH := $(FCM_ROOTDIR)
export FCM_CACHEPATH := $(FCM_CACHEDIR)
export FCM_CFGPATH := $(FCM_CFGDIR)
export FCM_SRCPATH := $(FCM_SRCDIR)
export FCM_BINPATH := $(FCM_BINDIR)
export FCM_ETCPATH := $(FCM_ETCDIR)
export FCM_LIBPATH := $(FCM_LIBDIR)
export FCM_TMPPATH := $(FCM_TMPDIR)
export FCM_DONEPATH := $(FCM_DONEDIR)
export FCM_FLAGSPATH := $(FCM_FLAGSDIR)
export FCM_INCPATH := $(FCM_INCDIR)
export FCM_PPSRCPATH := $(FCM_PPSRCDIR)
export FCM_OBJPATH := $(FCM_OBJDIR)

export PERL5LIB := /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/ext/FCM/bin/../lib

export AR = ar
export ARFLAGS = -r
export CC = cc
export CC_COMPILE = -c
export CC_DEFINE = -D
export CC_INCLUDE = -I
export CC_OUTPUT = -o
export CFLAGS = 
export CPPKEYS = 
export DIFF3FLAGS = -E -m
export FC = ftn
export FC_COMPILE = -c
export FC_DEFINE = -D
export FC_INCLUDE = -I
export FC_MODSEARCH = 
export FC_OUTPUT = -o
export FFLAGS = -s real64 -s integer32 -Ovector0 -hfp0 -O0 -hflex_mp=intolerant -e CID -I/opt/cray/netcdf-hdf5parallel/4.3.2/CRAY/83/include -I/opt/cray/hdf5-parallel/1.8.13/CRAY/83/include -I/common/moci/modules/packages/XIOS/2.0/24708/952/inc -I/common/moci/modules/packages/oasis3-mct/2.0/24708/1765/build/lib/mct -I/common/moci/modules/packages/oasis3-mct/2.0/24708/1765/build/lib/psmile.MPI1
export FPPKEYS = 
export LD = ftn
export LDFLAGS = -hbyteswapio -L/opt/cray/hdf5-parallel/1.8.13/CRAY/83/lib -L/opt/cray/netcdf-hdf5parallel/4.3.2/CRAY/83/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz -L/common/moci/modules/packages/XIOS/2.0/24708/952/lib  -L/common/moci/modules/packages/oasis3-mct/2.0/24708/1765/lib
export LD_LIBLINK = -l
export LD_LIBSEARCH = -L
export LD_OUTPUT = -o
export SHELL = /bin/sh

vpath % $(FCM_BINPATH)
vpath %.done $(FCM_DONEPATH)
vpath %.idone $(FCM_DONEPATH)
vpath %.etc $(FCM_ETCPATH)
vpath %.flags $(FCM_FLAGSPATH)
vpath % $(FCM_INCPATH)
vpath %.a $(FCM_LIBPATH)
vpath %.o $(FCM_OBJPATH)

FCM_BLD_TARGETS = make_domain_cfg.exe

.PHONY : all

all : $(FCM_BLD_TARGETS)

$(FCM_DONEDIR)/FCM_CP.dummy:
	touch $@

CFLAGS__nemo.flags: CFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

CPPKEYS__nemo.flags: CPPKEYS.flags
	touch $(FCM_FLAGSDIR)/$@

FFLAGS__nemo.flags: FFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

FPPKEYS__nemo.flags: FPPKEYS.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__nemo.flags: LDFLAGS.flags
	touch $(FCM_FLAGSDIR)/$@

LD__nemo.flags: LD.flags
	touch $(FCM_FLAGSDIR)/$@

OBJECTS__nemo = calendar.o closea.o daymod.o defprec.o depth_e3.o dom_oce.o domain.o domcfg.o domhgr.o dommsk.o domngb.o domstp.o domvvl.o domwri.o domzgr.o dtatsd.o eosbn2.o errioipsl.o flincom.o fliocom.o getincom.o histcom.o in_out_manager.o ioipsl.o iom.o iom_def.o iom_nf90.o iscplhsb.o iscplini.o iscplrst.o istate.o lbclnk.o lbcnfd.o lib_cray.o lib_fortran.o lib_mpp.o mathelp.o mes.o mppini.o nc4interface.o nemogcm.o oce.o par_kind.o par_oce.o phycst.o prtctl.o restart.o restcom.o step_oce.o stpctl.o stringop.o timing.o trc_oce.o wrk_nemo.o zgrloc.o zgrmes.o

libnemo.a: $(OBJECTS__nemo)
	fcm_internal archive $@ $^

CC.flags:
	touch $(FCM_FLAGSDIR)/$@

CFLAGS.flags: CC.flags
	touch $(FCM_FLAGSDIR)/$@

CPPKEYS.flags:
	touch $(FCM_FLAGSDIR)/$@

FC.flags:
	touch $(FCM_FLAGSDIR)/$@

FFLAGS.flags: FC.flags
	touch $(FCM_FLAGSDIR)/$@

FPPKEYS.flags:
	touch $(FCM_FLAGSDIR)/$@

LD.flags:
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS.flags:
	touch $(FCM_FLAGSDIR)/$@

export OBJECTS = $(OBJECTS__nemo)

libfcm_default.a: $(OBJECTS)
	fcm_internal archive $@ $^

FFLAGS__nemo__calendar.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

calendar.done: calendar.o errioipsl.done stringop.done
	touch $(FCM_DONEDIR)/$@

calendar.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/calendar.f90 FFLAGS__nemo__calendar.flags errioipsl.o stringop.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__closea.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

closea.done: closea.o dom_oce.done in_out_manager.done lbclnk.done lib_fortran.done lib_mpp.done oce.done phycst.done timing.done
	touch $(FCM_DONEDIR)/$@

closea.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/closea.f90 FFLAGS__nemo__closea.flags dom_oce.o in_out_manager.o lbclnk.o lib_fortran.o lib_mpp.o oce.o phycst.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__daymod.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

daymod.done: daymod.o dom_oce.done in_out_manager.done ioipsl.done iom.done phycst.done prtctl.done timing.done
	touch $(FCM_DONEDIR)/$@

daymod.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/daymod.f90 FFLAGS__nemo__daymod.flags dom_oce.o in_out_manager.o ioipsl.o iom.o phycst.o prtctl.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__defprec.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

defprec.done: defprec.o
	touch $(FCM_DONEDIR)/$@

defprec.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/defprec.f90 FFLAGS__nemo__defprec.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__depth_e3.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

FPPKEYS__nemo__depth_e3.flags: FPPKEYS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

depth_e3.done: depth_e3.o dom_oce.done in_out_manager.done lbclnk.done lib_mpp.done vectopt_loop_substitute.h90.idone
	touch $(FCM_DONEDIR)/$@

depth_e3.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/depth_e3.F90 FFLAGS__nemo__depth_e3.flags FPPKEYS__nemo__depth_e3.flags dom_oce.o in_out_manager.o lbclnk.o lib_mpp.o vectopt_loop_substitute.h90
	fcm_internal compile:F nemo $< $@ 1

FFLAGS__nemo__dom_oce.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dom_oce.done: dom_oce.o par_oce.done
	touch $(FCM_DONEDIR)/$@

dom_oce.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/dom_oce.f90 FFLAGS__nemo__dom_oce.flags par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__domain.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

domain.done: domain.o closea.done dom_oce.done domhgr.done dommsk.done domstp.done domvvl.done domwri.done domzgr.done in_out_manager.done ioipsl.done iom.done lbclnk.done lib_mpp.done oce.done phycst.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

domain.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/domain.f90 FFLAGS__nemo__domain.flags closea.o dom_oce.o domhgr.o dommsk.o domstp.o domvvl.o domwri.o domzgr.o in_out_manager.o ioipsl.o iom.o lbclnk.o lib_mpp.o oce.o phycst.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__domcfg.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

domcfg.done: domcfg.o dom_oce.done in_out_manager.done lib_mpp.done phycst.done timing.done
	touch $(FCM_DONEDIR)/$@

domcfg.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/domcfg.f90 FFLAGS__nemo__domcfg.flags dom_oce.o in_out_manager.o lib_mpp.o phycst.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__domhgr.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

domhgr.done: domhgr.o dom_oce.done domwri.done in_out_manager.done iom.done lib_mpp.done phycst.done timing.done
	touch $(FCM_DONEDIR)/$@

domhgr.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/domhgr.f90 FFLAGS__nemo__domhgr.flags dom_oce.o domwri.o in_out_manager.o iom.o lib_mpp.o phycst.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dommsk.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dommsk.done: dommsk.o dom_oce.done in_out_manager.done lbclnk.done lib_mpp.done oce.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

dommsk.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/dommsk.f90 FFLAGS__nemo__dommsk.flags dom_oce.o in_out_manager.o lbclnk.o lib_mpp.o oce.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__domngb.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

domngb.done: domngb.o dom_oce.done in_out_manager.done lib_mpp.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

domngb.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/domngb.f90 FFLAGS__nemo__domngb.flags dom_oce.o in_out_manager.o lib_mpp.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__domstp.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

domstp.done: domstp.o dom_oce.done in_out_manager.done lib_mpp.done oce.done
	touch $(FCM_DONEDIR)/$@

domstp.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/domstp.f90 FFLAGS__nemo__domstp.flags dom_oce.o in_out_manager.o lib_mpp.o oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__domvvl.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

domvvl.done: domvvl.o dom_oce.done in_out_manager.done iom.done lbclnk.done lib_mpp.done oce.done phycst.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

domvvl.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/domvvl.f90 FFLAGS__nemo__domvvl.flags dom_oce.o in_out_manager.o iom.o lbclnk.o lib_mpp.o oce.o phycst.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__domwri.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

domwri.done: domwri.o dom_oce.done in_out_manager.done iom.done lbclnk.done lib_mpp.done phycst.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

domwri.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/domwri.f90 FFLAGS__nemo__domwri.flags dom_oce.o in_out_manager.o iom.o lbclnk.o lib_mpp.o phycst.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__domzgr.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

domzgr.done: domzgr.o closea.done dom_oce.done in_out_manager.done iom.done lbclnk.done lib_mpp.done oce.done timing.done wrk_nemo.done zgrmes.done
	touch $(FCM_DONEDIR)/$@

domzgr.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/domzgr.f90 FFLAGS__nemo__domzgr.flags closea.o dom_oce.o in_out_manager.o iom.o lbclnk.o lib_mpp.o oce.o timing.o wrk_nemo.o zgrmes.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__dtatsd.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

dtatsd.done: dtatsd.o dom_oce.done fldread.done in_out_manager.done lib_mpp.done oce.done phycst.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

dtatsd.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/dtatsd.f90 FFLAGS__nemo__dtatsd.flags dom_oce.o fldread.o in_out_manager.o lib_mpp.o oce.o phycst.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__eosbn2.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

eosbn2.done: eosbn2.o dom_oce.done in_out_manager.done lbclnk.done lib_fortran.done lib_mpp.done phycst.done prtctl.done timing.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

eosbn2.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/eosbn2.f90 FFLAGS__nemo__eosbn2.flags dom_oce.o in_out_manager.o lbclnk.o lib_fortran.o lib_mpp.o phycst.o prtctl.o timing.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__errioipsl.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

errioipsl.done: errioipsl.o
	touch $(FCM_DONEDIR)/$@

errioipsl.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/errioipsl.f90 FFLAGS__nemo__errioipsl.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__flincom.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

flincom.done: flincom.o calendar.done errioipsl.done stringop.done
	touch $(FCM_DONEDIR)/$@

flincom.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/flincom.f90 FFLAGS__nemo__flincom.flags calendar.o errioipsl.o stringop.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__fliocom.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

fliocom.done: fliocom.o calendar.done defprec.done errioipsl.done stringop.done
	touch $(FCM_DONEDIR)/$@

fliocom.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/fliocom.f90 FFLAGS__nemo__fliocom.flags calendar.o defprec.o errioipsl.o stringop.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__getincom.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

getincom.done: getincom.o errioipsl.done stringop.done
	touch $(FCM_DONEDIR)/$@

getincom.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/getincom.f90 FFLAGS__nemo__getincom.flags errioipsl.o stringop.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__histcom.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

histcom.done: histcom.o calendar.done errioipsl.done fliocom.done mathelp.done nc4interface.done stringop.done
	touch $(FCM_DONEDIR)/$@

histcom.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/histcom.f90 FFLAGS__nemo__histcom.flags calendar.o errioipsl.o fliocom.o mathelp.o nc4interface.o stringop.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__in_out_manager.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

in_out_manager.done: in_out_manager.o nc4interface.done par_oce.done
	touch $(FCM_DONEDIR)/$@

in_out_manager.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/in_out_manager.f90 FFLAGS__nemo__in_out_manager.flags nc4interface.o par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__ioipsl.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

ioipsl.done: ioipsl.o calendar.done errioipsl.done flincom.done fliocom.done getincom.done histcom.done mathelp.done restcom.done stringop.done
	touch $(FCM_DONEDIR)/$@

ioipsl.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/ioipsl.f90 FFLAGS__nemo__ioipsl.flags calendar.o errioipsl.o flincom.o fliocom.o getincom.o histcom.o mathelp.o restcom.o stringop.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__iom.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

iom.done: iom.o dom_oce.done in_out_manager.done iom_def.done iom_nf90.done lbclnk.done lib_mpp.done
	touch $(FCM_DONEDIR)/$@

iom.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/iom.f90 FFLAGS__nemo__iom.flags dom_oce.o in_out_manager.o iom_def.o iom_nf90.o lbclnk.o lib_mpp.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__iom_def.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

iom_def.done: iom_def.o par_kind.done
	touch $(FCM_DONEDIR)/$@

iom_def.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/iom_def.f90 FFLAGS__nemo__iom_def.flags par_kind.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__iom_nf90.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

iom_nf90.done: iom_nf90.o dom_oce.done in_out_manager.done iom_def.done lbclnk.done lib_mpp.done
	touch $(FCM_DONEDIR)/$@

iom_nf90.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/iom_nf90.f90 FFLAGS__nemo__iom_nf90.flags dom_oce.o in_out_manager.o iom_def.o lbclnk.o lib_mpp.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__iscplhsb.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

iscplhsb.done: iscplhsb.o dom_oce.done domngb.done domwri.done in_out_manager.done iscplini.done lbclnk.done lib_fortran.done lib_mpp.done oce.done phycst.done sbc_oce.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

iscplhsb.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/iscplhsb.f90 FFLAGS__nemo__iscplhsb.flags dom_oce.o domngb.o domwri.o in_out_manager.o iscplini.o lbclnk.o lib_fortran.o lib_mpp.o oce.o phycst.o sbc_oce.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__iscplini.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

iscplini.done: iscplini.o dom_oce.done in_out_manager.done lib_fortran.done lib_mpp.done oce.done
	touch $(FCM_DONEDIR)/$@

iscplini.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/iscplini.f90 FFLAGS__nemo__iscplini.flags dom_oce.o in_out_manager.o lib_fortran.o lib_mpp.o oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__iscplrst.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

iscplrst.done: iscplrst.o dom_oce.done domvvl.done domwri.done in_out_manager.done iom.done iscplhsb.done iscplini.done lbclnk.done lib_fortran.done lib_mpp.done oce.done phycst.done sbc_oce.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

iscplrst.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/iscplrst.f90 FFLAGS__nemo__iscplrst.flags dom_oce.o domvvl.o domwri.o in_out_manager.o iom.o iscplhsb.o iscplini.o lbclnk.o lib_fortran.o lib_mpp.o oce.o phycst.o sbc_oce.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__istate.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

istate.done: istate.o c1d.done daymod.done divhor.done dom_oce.done domvvl.done dtatsd.done dtauvd.done eosbn2.done in_out_manager.done iom.done iscplrst.done lbclnk.done ldftra.done lib_mpp.done oce.done phycst.done restart.done timing.done wrk_nemo.done zdf_oce.done
	touch $(FCM_DONEDIR)/$@

istate.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/istate.f90 FFLAGS__nemo__istate.flags c1d.o daymod.o divhor.o dom_oce.o domvvl.o dtatsd.o dtauvd.o eosbn2.o in_out_manager.o iom.o iscplrst.o lbclnk.o ldftra.o lib_mpp.o oce.o phycst.o restart.o timing.o wrk_nemo.o zdf_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__lbclnk.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

lbclnk.done: lbclnk.o lib_mpp.done
	touch $(FCM_DONEDIR)/$@

lbclnk.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/lbclnk.f90 FFLAGS__nemo__lbclnk.flags lib_mpp.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__lbcnfd.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

lbcnfd.done: lbcnfd.o dom_oce.done in_out_manager.done
	touch $(FCM_DONEDIR)/$@

lbcnfd.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/lbcnfd.f90 FFLAGS__nemo__lbcnfd.flags dom_oce.o in_out_manager.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__lib_cray.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

lib_cray.done: lib_cray.o
	touch $(FCM_DONEDIR)/$@

lib_cray.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/lib_cray.f90 FFLAGS__nemo__lib_cray.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__lib_fortran.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

lib_fortran.done: lib_fortran.o dom_oce.done in_out_manager.done lib_mpp.done par_oce.done
	touch $(FCM_DONEDIR)/$@

lib_fortran.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/lib_fortran.f90 FFLAGS__nemo__lib_fortran.flags dom_oce.o in_out_manager.o lib_mpp.o par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__lib_mpp.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

lib_mpp.done: lib_mpp.o dom_oce.done in_out_manager.done lbcnfd.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

lib_mpp.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/lib_mpp.f90 FFLAGS__nemo__lib_mpp.flags dom_oce.o in_out_manager.o lbcnfd.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__make_domain_cfg.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

LDFLAGS__nemo__make_domain_cfg.flags: LDFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

LD__nemo__make_domain_cfg.flags: LD__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

make_domain_cfg.exe: make_domain_cfg.o LD__nemo__make_domain_cfg.flags LDFLAGS__nemo__make_domain_cfg.flags nemogcm.done
	fcm_internal load:F nemo $< $@

make_domain_cfg.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/make_domain_cfg.f90 FFLAGS__nemo__make_domain_cfg.flags nemogcm.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__mathelp.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

mathelp.done: mathelp.o errioipsl.done stringop.done
	touch $(FCM_DONEDIR)/$@

mathelp.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/mathelp.f90 FFLAGS__nemo__mathelp.flags errioipsl.o stringop.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__mes.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

FPPKEYS__nemo__mes.flags: FPPKEYS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

mes.done: mes.o closea.done dom_oce.done in_out_manager.done iom.done lbclnk.done lib_mpp.done oce.done timing.done vectopt_loop_substitute.h90.idone wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

mes.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/mes.F90 FFLAGS__nemo__mes.flags FPPKEYS__nemo__mes.flags closea.o dom_oce.o in_out_manager.o iom.o lbclnk.o lib_mpp.o oce.o timing.o vectopt_loop_substitute.h90 wrk_nemo.o
	fcm_internal compile:F nemo $< $@ 1

FFLAGS__nemo__mppini.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

mppini.done: mppini.o dom_oce.done in_out_manager.done ioipsl.done iom.done lib_mpp.done
	touch $(FCM_DONEDIR)/$@

mppini.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/mppini.f90 FFLAGS__nemo__mppini.flags dom_oce.o in_out_manager.o ioipsl.o iom.o lib_mpp.o
	fcm_internal compile:F nemo $< $@

mppini_2.h90: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/mppini_2.h90
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

mppini_2.h90.idone: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/mppini_2.h90
	touch $(FCM_DONEDIR)/$@

FFLAGS__nemo__nc4interface.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

nc4interface.done: nc4interface.o
	touch $(FCM_DONEDIR)/$@

nc4interface.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/nc4interface.f90 FFLAGS__nemo__nc4interface.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__nemogcm.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

nemogcm.done: nemogcm.o dom_oce.done domain.done domcfg.done lbcnfd.done lib_fortran.done lib_mpp.done mppini.done phycst.done step_oce.done
	touch $(FCM_DONEDIR)/$@

nemogcm.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/nemogcm.f90 FFLAGS__nemo__nemogcm.flags dom_oce.o domain.o domcfg.o lbcnfd.o lib_fortran.o lib_mpp.o mppini.o phycst.o step_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__oce.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

oce.done: oce.o lib_mpp.done par_oce.done
	touch $(FCM_DONEDIR)/$@

oce.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/oce.f90 FFLAGS__nemo__oce.flags lib_mpp.o par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__par_kind.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

par_kind.done: par_kind.o
	touch $(FCM_DONEDIR)/$@

par_kind.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/par_kind.f90 FFLAGS__nemo__par_kind.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__par_oce.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

par_oce.done: par_oce.o par_kind.done
	touch $(FCM_DONEDIR)/$@

par_oce.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/par_oce.f90 FFLAGS__nemo__par_oce.flags par_kind.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__phycst.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

phycst.done: phycst.o in_out_manager.done par_oce.done
	touch $(FCM_DONEDIR)/$@

phycst.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/phycst.f90 FFLAGS__nemo__phycst.flags in_out_manager.o par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__prtctl.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

prtctl.done: prtctl.o dom_oce.done in_out_manager.done lib_mpp.done wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

prtctl.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/prtctl.f90 FFLAGS__nemo__prtctl.flags dom_oce.o in_out_manager.o lib_mpp.o wrk_nemo.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__restart.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

restart.done: restart.o diurnal_bulk.done dom_oce.done eosbn2.done in_out_manager.done iom.done oce.done phycst.done sbc_ice.done trdmxl_oce.done
	touch $(FCM_DONEDIR)/$@

restart.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/restart.f90 FFLAGS__nemo__restart.flags diurnal_bulk.o dom_oce.o eosbn2.o in_out_manager.o iom.o oce.o phycst.o sbc_ice.o trdmxl_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__restcom.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

restcom.done: restcom.o calendar.done errioipsl.done fliocom.done mathelp.done stringop.done
	touch $(FCM_DONEDIR)/$@

restcom.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/restcom.f90 FFLAGS__nemo__restcom.flags calendar.o errioipsl.o fliocom.o mathelp.o stringop.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__step_oce.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

step_oce.done: step_oce.o daymod.done dom_oce.done eosbn2.done in_out_manager.done iom.done lbclnk.done oce.done prtctl.done timing.done
	touch $(FCM_DONEDIR)/$@

step_oce.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/step_oce.f90 FFLAGS__nemo__step_oce.flags daymod.o dom_oce.o eosbn2.o in_out_manager.o iom.o lbclnk.o oce.o prtctl.o timing.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__stpctl.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

stpctl.done: stpctl.o c1d.done dom_oce.done in_out_manager.done lbclnk.done lib_mpp.done oce.done
	touch $(FCM_DONEDIR)/$@

stpctl.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/stpctl.f90 FFLAGS__nemo__stpctl.flags c1d.o dom_oce.o in_out_manager.o lbclnk.o lib_mpp.o oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__stringop.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

stringop.done: stringop.o
	touch $(FCM_DONEDIR)/$@

stringop.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/stringop.f90 FFLAGS__nemo__stringop.flags
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__timing.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

timing.done: timing.o dom_oce.done in_out_manager.done lib_mpp.done
	touch $(FCM_DONEDIR)/$@

timing.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/timing.f90 FFLAGS__nemo__timing.flags dom_oce.o in_out_manager.o lib_mpp.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__trc_oce.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

trc_oce.done: trc_oce.o dom_oce.done in_out_manager.done lib_mpp.done par_oce.done
	touch $(FCM_DONEDIR)/$@

trc_oce.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/trc_oce.f90 FFLAGS__nemo__trc_oce.flags dom_oce.o in_out_manager.o lib_mpp.o par_oce.o
	fcm_internal compile:F nemo $< $@

vectopt_loop_substitute.h90: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/vectopt_loop_substitute.h90
	cp $< $(FCM_INCDIR)
	chmod u+w $(FCM_INCDIR)/$@

vectopt_loop_substitute.h90.idone: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/vectopt_loop_substitute.h90
	touch $(FCM_DONEDIR)/$@

FFLAGS__nemo__wrk_nemo.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

wrk_nemo.done: wrk_nemo.o par_oce.done
	touch $(FCM_DONEDIR)/$@

wrk_nemo.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/wrk_nemo.f90 FFLAGS__nemo__wrk_nemo.flags par_oce.o
	fcm_internal compile:F nemo $< $@

FFLAGS__nemo__zgrloc.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

FPPKEYS__nemo__zgrloc.flags: FPPKEYS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

zgrloc.done: zgrloc.o closea.done depth_e3.done dom_oce.done in_out_manager.done iom.done lbclnk.done lib_mpp.done mes.done oce.done timing.done vectopt_loop_substitute.h90.idone wrk_nemo.done
	touch $(FCM_DONEDIR)/$@

zgrloc.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/zgrloc.F90 FFLAGS__nemo__zgrloc.flags FPPKEYS__nemo__zgrloc.flags closea.o depth_e3.o dom_oce.o in_out_manager.o iom.o lbclnk.o lib_mpp.o mes.o oce.o timing.o vectopt_loop_substitute.h90 wrk_nemo.o
	fcm_internal compile:F nemo $< $@ 1

FFLAGS__nemo__zgrmes.flags: FFLAGS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

FPPKEYS__nemo__zgrmes.flags: FPPKEYS__nemo.flags
	touch $(FCM_FLAGSDIR)/$@

zgrmes.done: zgrmes.o closea.done depth_e3.done dom_oce.done in_out_manager.done iom.done lbclnk.done lib_mpp.done mes.done oce.done timing.done vectopt_loop_substitute.h90.idone wrk_nemo.done zgrloc.done
	touch $(FCM_DONEDIR)/$@

zgrmes.o: /projects/jmmp/dbruciaf/NEMO/CHECKOUTS_RELEASE-4.0-HEAD/NEMO_4.0-HEAD_test_MEs/tools/DOMAINcfg/src/zgrmes.F90 FFLAGS__nemo__zgrmes.flags FPPKEYS__nemo__zgrmes.flags closea.o depth_e3.o dom_oce.o in_out_manager.o iom.o lbclnk.o lib_mpp.o mes.o oce.o timing.o vectopt_loop_substitute.h90 wrk_nemo.o zgrloc.o
	fcm_internal compile:F nemo $< $@ 1

# EOF
