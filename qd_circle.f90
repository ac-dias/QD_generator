program main

	implicit none

	double precision :: a
	double precision,parameter :: theta=0.710

	double precision,dimension(2) :: r1,r2,rp,r
	integer :: m,n,i,j,erro

	double precision :: rcut !tamanho do ponto quantico
	double precision :: modv
	double precision,parameter :: tol=1.5
	
	integer :: natom

	character(len=3) :: tm,ch
	character(len=50) :: outfile


	write(*,*) "Lattice Constant in Angstrom"
	read(*,*) a
	
	write(*,*) "Quantum Dot Radius in Angstrom"
	read(*,*) rcut
	
	write(outfile,"(a3,I0,a4)") 'qd_',int(rcut),".xyz"
	
	OPEN(UNIT=300, FILE= outfile,STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de saida xyz"		
	
	write(*,*) "Transition metal species"
	read(*,*) tm

	write(*,*) "chalcogen/Halogen species"
	read(*,*) ch



	r1(1)= a
	r1(2)= 0.

	r2(1)= a/2.
	r2(2)= (sqrt(3.)/2.)*a

	rp(1)= 0.
	rp(2)= a/sqrt(3.)
	
	call atomcount(a,rcut,natom)
	
	write(300,*) natom
	write(300,*) "qdot"

	!call circle(rcut)
	!construindo a rede m

	do m=-100,100

		do n=-100,100

			r(1)= dble(m)*r1(1)+dble(n)*r2(1)
			r(2)=dble(m)*r1(2)+dble(n)*r2(2)


			call modvec(r,modv)


			if (modv .lt.  rcut+tol ) then
	

				write(300,*) trim(tm),"      ",r(1),r(2),0.0

			else 

				continue

			end if
			

		end do


	end do

	!construindo a rede x

	do m=-100,100

		do n=-100,100

			r(1)= dble(m)*r1(1)+dble(n)*r2(1)+rp(1)
			r(2)=dble(m)*r1(2)+dble(n)*r2(2)+rp(2)


			call modvec(r,modv)


			if (modv .lt.  rcut+tol ) then
	

				write(300,*) trim(ch),"      ",r(1),r(2), (a/(sqrt(3.)*cos(theta)))*sin(theta)
				write(300,*) trim(ch),"      ",r(1),r(2),-(a/(sqrt(3.)*cos(theta)))*sin(theta)


			else 

				continue

			end if
			

		end do


	end do
	
	close(300)


end program main

subroutine atomcount(a,rcut,natom)

	implicit none
	
	double precision :: a
	double precision,parameter :: theta=0.710

	double precision,dimension(2) :: r1,r2,rp,r
	integer :: m,n,i,j,erro

	double precision :: rcut !tamanho do ponto quantico
	double precision :: modv
	double precision,parameter :: tol=1.5
	
	integer :: natom
	
	r1(1)= a
	r1(2)= 0.

	r2(1)= a/2.
	r2(2)= (sqrt(3.)/2.)*a

	rp(1)= 0.
	rp(2)= a/sqrt(3.)	
	
	natom = 0
	
	do m=-100,100

		do n=-100,100

			r(1)= dble(m)*r1(1)+dble(n)*r2(1)
			r(2)=dble(m)*r1(2)+dble(n)*r2(2)


			call modvec(r,modv)


			if (modv .lt.  rcut+tol ) then
	

				natom = natom + 1

			else 

				continue

			end if
			

		end do


	end do

	!construindo a rede x

	do m=-100,100

		do n=-100,100

			r(1)= dble(m)*r1(1)+dble(n)*r2(1)+rp(1)
			r(2)=dble(m)*r1(2)+dble(n)*r2(2)+rp(2)


			call modvec(r,modv)


			if (modv .lt.  rcut+tol ) then
	

				natom = natom + 2


			else 

				continue

			end if
			

		end do


	end do	
	

end subroutine atomcount

subroutine modvec(r,modv)

	implicit none

	double precision,dimension(2) :: r
	double precision :: modv

	modv=sqrt(r(1)*r(1)+r(2)*r(2))


end subroutine modvec


