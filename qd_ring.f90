program main

	implicit none

	real :: a
	real,parameter :: theta=0.710

	real,dimension(2) :: r1,r2,rp,r
	integer :: m,n,i,j,erro,natom

	real :: rcutex !tamanho do ponto quantico
	real :: rcutin !tamanho do ponto quantico
	real :: modv
	
	character(len=3) :: tm,ch
	character(len=50) :: outfile	
   	

    	
	write(*,*) "Lattice Constant in Angstrom"
	read(*,*) a
	
	write(*,*) "Quantum Dot External Radius in Angstrom"
	read(*,*) rcutex
	
	write(*,*) "Quantum Dot Internal Radius in Angstrom"
	read(*,*) rcutin	
	
    	if (rcutin .ge. rcutex) then
    	
    	 stop "Internal Radius higher than External Radius"	
    	
    	end if	
	
	write(outfile,"(a6,I0,a1,I0,a4)") 'qring_',int(rcutex),"_",int(rcutin),".xyz"
	
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

	call atomcount(a,rcutex,rcutin,natom)
	
	write(300,*) natom
	write(300,*) "qdot"

	!construindo a rede m

	do m=-100,100

		do n=-100,100

			r(1)= m*r1(1)+n*r2(1)
			r(2)=m*r1(2)+n*r2(2)


			call modvec(r,modv)


			if ( (modv .le. rcutex) .and. (modv .ge. rcutin)) then
	
				write(300,*) "Mo","      ",r(1),r(2),0.0

			else 

				continue

			end if
			

		end do


	end do

	!construindo a rede x

	do m=-100,100

		do n=-100,100

			r(1)= m*r1(1)+n*r2(1)+rp(1)
			r(2)=m*r1(2)+n*r2(2)+rp(2)


			call modvec(r,modv)


			if ( (modv .le. rcutex) .and. (modv .ge. rcutin)) then
	
				write(300,*) "S","      ",r(1),r(2), (a/(sqrt(3.)*cos(theta)))*sin(theta)
				write(300,*) "S","      ",r(1),r(2),-(a/(sqrt(3.)*cos(theta)))*sin(theta)

			else 

				continue

			end if
			

		end do


	end do


	close(300)



end program main

subroutine atomcount(a,rcutex,rcutin,natom)

	implicit none
	
	real :: a
	real,parameter :: theta=0.710

	real,dimension(2) :: r1,r2,rp,r
	integer :: m,n,i,j,erro

	real :: rcutex !tamanho do ponto quantico
	real :: rcutin !tamanho do ponto quantico
	real :: modv
	real,parameter :: tol=1.5
	
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

			r(1)= m*r1(1)+n*r2(1)
			r(2)=m*r1(2)+n*r2(2)


			call modvec(r,modv)


			if ( (modv .le. rcutex) .and. (modv .ge. rcutin)) then
	
				natom = natom +1

			else 

				continue

			end if
			

		end do


	end do

	!construindo a rede x

	do m=-100,100

		do n=-100,100

			r(1)= m*r1(1)+n*r2(1)+rp(1)
			r(2)=m*r1(2)+n*r2(2)+rp(2)


			call modvec(r,modv)


			if ( (modv .le. rcutex) .and. (modv .ge. rcutin)) then
	
				natom = natom +2

			else 

				continue

			end if
			

		end do


	end do
	

end subroutine atomcount

subroutine modvec(r,modv)

	implicit none

	real,dimension(2) :: r
	real :: modv

	modv=sqrt(r(1)*r(1)+r(2)*r(2))


end subroutine modvec

subroutine ring(rcutin,rcutex)


	implicit none

	real :: rcutin,rcutex
	real, parameter :: pi=acos(-1.)
	real :: phi
	integer,parameter :: jtot=1000

	integer :: erro,j

	OPEN(UNIT=300, FILE='ring.dat',STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de saida da rede m"

	do j=1,jtot

		phi=(2.*pi)*((j-1.)/(jtot-1.))

		write(300,*) rcutin*cos(phi),rcutin*sin(phi)

	end do

	do j=1,jtot

		phi=(2.*pi)*((j-1.)/(jtot-1.))

		write(300,*) rcutex*cos(phi),rcutex*sin(phi)

	end do


	close(300)

end subroutine ring
