include ../standard_defs.mk

all: mod_netcdf_error.o *.o move

mod_netcdf_error.o : mod_netcdf_error.f
	$(CC) -c $(LIBS) $(INCLUDES) mod_netcdf_error.f $(LINKS)

*.o : *.f
	$(CC) -c $(LIBS) $(INCLUDES) *.f $(LINKS)
move:
	mv *.o *.mod $(BUILD_DIR)

clean:
	rm *.o *.mod 	

