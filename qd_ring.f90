program main

	implicit none

	real,parameter :: a=3.166
	real,parameter :: theta=0.710

	real,dimension(2) :: r1,r2,rp,r
	integer :: m,n,i,j,erro

	real,parameter :: rcutex=50. !tamanho do ponto quantico
	real,parameter :: rcutin=20. !tamanho do ponto quantico
	real :: modv

	OPEN(UNIT=200, FILE='qdring2.xyz',STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de saida da rede m"


	r1(1)= a
	r1(2)= 0.

	r2(1)= a/2.
	r2(2)= (sqrt(3.)/2.)*a

	rp(1)= 0.
	rp(2)= a/sqrt(3.)

	call ring(rcutin,rcutex) !desenha as bordas do anel

	!construindo a rede m

	do m=-100,100

		do n=-100,100

			r(1)= m*r1(1)+n*r2(1)
			r(2)=m*r1(2)+n*r2(2)


			call modvec(r,modv)


			if ( (modv .le. rcutex) .and. (modv .ge. rcutin)) then
	
				write(200,*) "Mo","      ",r(1),r(2),0.0

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
	
				write(200,*) "S","      ",r(1),r(2), (a/(sqrt(3.)*cos(theta)))*sin(theta)
				write(200,*) "S","      ",r(1),r(2),-(a/(sqrt(3.)*cos(theta)))*sin(theta)

			else 

				continue

			end if
			

		end do


	end do


	close(200)
	close(201)


end program main

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
