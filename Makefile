CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall -O2
CXXFLAGS=	$(CFLAGS)
CPPFLAGS=
INCLUDES=
OBJS=		lv89.o edlib.o
PROG=		ed-test
LIBS=		-lz -lpthread -lm
LIBS_WFA2=

ifneq ($(WFA2_ROOT),)
	CPPFLAGS+=-D_USE_WFA2
	LIBS_WFA2=-L$(WFA2_ROOT)/lib -lwfa
	INCLUDES+=-I$(WFA2_ROOT)
endif

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address
endif

.SUFFIXES:.c .cpp .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

$(PROG):$(OBJS) main.o
		$(CXX) $(CFLAGS) $^ -o $@ $(LIBS_WFA2) $(LIBS)

clean:
		rm -fr gmon.out *.o a.out $(PROG) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c *.cpp)

# DO NOT DELETE

example.o: lv89.h
lv89-full.o: lv89.h
lv89-semi.o: lv89.h
lv89.o: lv89.h
main.o: lv89.h edlib.h ketopt.h kseq.h
edlib.o: edlib.h
