#------------------------------------------------------------------
OPT += -DNFW
#OPT += -DBURKERT                            ## for test-only
#------------------------------------------------------------------

CC       =   gcc   
OPTIMIZE =   -O3 -Wall 
GSL_INCL =  
GSL_LIBS =  

OPTIONS =  $(OPTIMIZE) $(OPT)

EXEC   = exe

OBJS   = main.o allvars.o init.o run.o mass_integral.o \
         finalize.o spline_alloc.o

INCL   = define.h allvars.h proto.h Makefile


CFLAGS = $(OPTIONS) $(GSL_INCL) 


LIBS   = $(GSL_LIBS) -lgsl -lgslcblas -lm 

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC) *.gch
