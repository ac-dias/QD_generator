include makefile.inc

all:

	$(FOR) ./qd_circle.f90 -o $(DIR)qd_circle.x 
	$(FOR) ./qd_hexagonal.f90 -o $(DIR)qd_hexagonal.x
	$(FOR) ./qd_rectangular.f90 -o $(DIR)qd_rectangular.x
	$(FOR) ./qd_ring.f90 -o $(DIR)qd_ring.x
	$(FOR) ./qd_tri-ac.f90 -o $(DIR)qd_tri-ac.x
	$(FOR) ./qd_tri-zz.f90 -o $(DIR)qd_tri-zz.x					

clean:
	@echo "--- Cleaning build, bin, and .mod files ---"
	@rm -rf $(DIR)/*.x




.PHONY: all clean 
 


