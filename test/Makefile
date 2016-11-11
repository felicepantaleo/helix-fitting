# make sure, env variable LD_LIBRARY_PATH = /usr/local/lib

#ML=$(shell [ -e Makefile.locale ]  && echo 1)
#ifeq ($(ML),1)
#	include Makefile.locale
#endif
	
CXXFLAGS += -DETCDIR=\"$(DESTDIR)/etc/\" $(SOLFLAG) $(INCLUDES) $(SETOPTIMIZE) -DCOLORS $(SETDEBUG) $(SETVERBOSE) -I/usr/local/include -I$(CLHEP)include/

OBJS = traxx.o
GLOBJS = gltraxx.o
CFILES = traxx.cpp

ifneq ($(DEBUG),0)
	override CXXFLAGS += -pg -g
endif

INCLUDES = -I/usr/X11R6/include $(HEPVIS_INCLUDE) \
   -I$(CLHEP)/include\
   -I../interface -I../hepvis/

gl : override CXXFLAGS += -fPIC -I$(INVENTOR)include -I/usr/local/include -DHEPVIS # -DGL
gl : override BINARY = gltraxx

nogl: $(OBJS)

nogl-nodebug: $(OBJS)
	
gltraxx.o: traxx.cpp
	$(CXX) $(CXXFLAGS) -c -o gltraxx.o traxx.cpp
gl: $(GLOBJS)
	
clean:
	rm -f $(OBJS) $(GLOBJS) *~

mrproper: clean
	rm -f .\#* .*.swp core .ignore_starter.m

ctags:
	ctags --no-warn --c++ *.cc *.h
	echo "!_TAG_FILE_SORTED	0" > /tmp/tags
	cat tags |  sed -e 's/^\w*:://' >> /tmp/tags
	cp /tmp/tags .

FORCE:
