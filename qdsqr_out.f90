program main

	implicit none

	double precision :: a
	double precision,parameter :: theta=0.710

	double precision,dimension(2) :: r1,r2,rp,r
	integer :: m,n,i,j,erro

	integer,parameter :: mmax=17 !tamanho do ponto quantico 
	double precision,parameter :: tol=0.1
	double precision :: modv

	character(*),parameter :: auxiliar='/home/adias/codigosalex/qdmx2/parametros-tmd/aux_mos2t.dat'

	character(*),parameter :: redem='rede_m_sqr.dat'
	character(*),parameter :: redex='rede_x_sqr.dat'
	character(*),parameter :: xyz='qdot_sqr.xyz'

	double precision :: alpha

	!input
	OPEN(UNIT=100, FILE= auxiliar,STATUS='old', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de saida da rede m"

	!output
	OPEN(UNIT=200, FILE= redem ,STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de saida da rede m"
	OPEN(UNIT=201, FILE= redex ,STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de saida da rede x"
	OPEN(UNIT=300, FILE= xyz,STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de saida xyz"

	read(100,*) a

	r1(1)= a
	r1(2)= 0.

	r2(1)= a/2.
	r2(2)= (sqrt(3.)/2.)*a

	rp(1)= 0.
	rp(2)= a/sqrt(3.)

	!call circle(rcut)
	!construindo a rede m

	do n=0,3000

		do m=-n,mmax/2


			

			r(1)= dble(2*m)*r1(1)+dble(n)*r2(1)
			r(2)=dble(2*m)*r1(2)+dble(n)*r2(2)


				!if ( (r(2) .gt.(a*mmax)) .or. (r(1) .gt.(a*mmax)) .or. (r(1) .lt. 0d0 ) ) go to 100
				if ( (r(2) .gt.(mmax*r2(2)+tol)) .or. (r(1) .gt.(a*mmax+tol)) .or. (r(1) .lt. 0d0 ) ) go to 100


	
				write(200,*) 2*m,n

				write(300,*) "Mo","      ",r(1),r(2),0.0

100				continue								
				

			r(1)= dble(2*m+1)*r1(1)+dble(n)*r2(1)
			r(2)=dble(2*m+1)*r1(2)+dble(n)*r2(2)

				!if ( (r(2) .gt.(a*mmax)) .or. (r(1) .gt.(a*mmax)) .or. (r(1) .lt. 0d0 ) ) go to 200
				if ( (r(2) .gt.(mmax*r2(2)+tol)) .or. (r(1) .gt.(a*mmax+tol)) .or. (r(1) .lt. 0d0 ) ) go to 200

	
				write(200,*) 2*m+1,n

				write(300,*) "Mo","      ",r(1),r(2),0.0


200				continue
			

		end do


	end do

	!construindo a rede x

	do n=-1,3000

		do m=-n-1,mmax/2

			!if ( (mmax .eq. 7) .or. (mmax .eq. 8)) then

			!	alpha=0d0

			!else 

			!	alpha=1d0

			!end if
			

			r(1)= dble(2*m)*r1(1)+dble(n)*r2(1)+rp(1)
			r(2)=dble(2*m)*r1(2)+dble(n)*r2(2)+rp(2)


				!if ( (r(2) .gt.(a*(mmax)-alpha*rp(2))) .or. (r(1) .gt.(a*mmax)) .or. (r(1) .lt. 0d0 ) )  go to 300
				if ( (r(2) .gt.(mmax*r2(2)+tol)) .or. (r(1) .gt.(a*mmax+tol)) .or. (r(1) .lt. 0d0 ) ) go to 300

	
				write(201,*) 2*m,n

				write(300,*) "S","      ",r(1),r(2), (a/(sqrt(3.)*cos(theta)))*sin(theta)
				write(300,*) "S","      ",r(1),r(2),-(a/(sqrt(3.)*cos(theta)))*sin(theta)

300				continue								
				

			r(1)= dble(2*m+1)*r1(1)+dble(n)*r2(1)+rp(1)
			r(2)=dble(2*m+1)*r1(2)+dble(n)*r2(2)+rp(2)

				!if ( (r(2) .gt.(a*(mmax)-alpha*rp(2))) .or. (r(1) .gt.(a*mmax)) .or. (r(1) .lt. 0d0 ) )  go to 400
				if ( (r(2) .gt.(mmax*r2(2)+tol)) .or. (r(1) .gt.(a*mmax+tol)) .or. (r(1) .lt. 0d0 ) ) go to 400
	
				write(201,*) 2*m+1,n

				write(300,*) "S","      ",r(1),r(2), (a/(sqrt(3.)*cos(theta)))*sin(theta)
				write(300,*) "S","      ",r(1),r(2),-(a/(sqrt(3.)*cos(theta)))*sin(theta)

400				continue
			

		end do


	end do




	close(100)
	close(200)
	close(201)
	close(300)


end program main

subroutine modvec(r,modv)

	implicit none

	double precision,dimension(2) :: r
	double precision :: modv

	modv=sqrt(r(1)*r(1)+r(2)*r(2))


end subroutine modvec

subroutine circle(rcut)


	implicit none

	double precision :: rcut
	double precision, parameter :: pi=acos(-1.)
	double precision :: phi
	integer,parameter :: jtot=1000

	integer :: erro,j

	OPEN(UNIT=300, FILE='circle.dat',STATUS='unknown', IOSTAT=erro)
    	if (erro/=0) stop "Erro na abertura do arquivo de saida da rede m"

	do j=1,jtot

		phi=(2.*pi)*((j-1.)/(jtot-1.))

		write(300,*) rcut*cos(phi),rcut*sin(phi)

	end do


	close(300)

end subroutine circle
