program main

	implicit none

	double precision :: a
	double precision,parameter :: theta=0.710

	double precision,dimension(2) :: r1,r2,rp,r
	integer :: m,n,i,j,erro

	integer,parameter :: mmax=6 !tamanho do ponto quantico
	double precision :: modv

	character(*),parameter :: auxiliar='/home/adias/codigosalex/qdmx2/parametros-tmd/aux_mos2t.dat'

	character(*),parameter :: redem='rede_m_tri-ac.dat'
	character(*),parameter :: redex='rede_x_tri-ac.dat'
	character(*),parameter :: xyz='qdot_tri-ac.xyz'

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

	do n=-1,mmax+1

		do m=-1,mmax-n+1

			if (n==-1 .and. m==-1)  cycle
			if (n==-1 .and. m==mmax-n+1) cycle

			r(1)= dble(m)*r1(1)+dble(n)*r2(1)
			r(2)=dble(m)*r1(2)+dble(n)*r2(2)


	
				write(201,*) m,n

				write(300,*) "S","      ",r(1),r(2), (a/(sqrt(3.)*cos(theta)))*sin(theta)
				write(300,*) "S","      ",r(1),r(2),-(a/(sqrt(3.)*cos(theta)))*sin(theta)

			

		end do


	end do

	!construindo a rede x

	do n=-1,mmax

		do m=0,mmax-n

			r(1)= dble(m)*r1(1)+dble(n)*r2(1)+rp(1)
			r(2)=dble(m)*r1(2)+dble(n)*r2(2)+rp(2)


	
				write(200,*) m,n

				write(300,*) "Mo","      ",r(1),r(2),0.0
		


			

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
