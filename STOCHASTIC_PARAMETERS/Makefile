include ../standard_defs.mk

all: mod_stochastic_parameters.o move

mod_stochastic_parameters.o : mod_stochastic_parameters.f
	$(CC) -c $(LIBS) $(INCLUDES) *.f $(LINKS)

move:
	mv *.o *.mod $(BUILD_DIR)

clean:
	rm *.o *.mod 	

