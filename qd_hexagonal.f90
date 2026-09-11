program main

	implicit none

	double precision :: a
	double precision,parameter :: theta=0.710

	double precision,dimension(2) :: r1,r2,rp,r
	integer :: m,n,i,j,erro,natom

	integer :: mmax !tamanho do ponto quantico
	double precision :: modv

	character(len=3) :: tm,ch
	character(len=50) :: outfile
	
	
	
	write(*,*) "Lattice Constant in Angstrom"
	read(*,*) a	
	
	write(*,*) "Quantum Size Factor (integer value)"
	read(*,*) mmax
			

	write(outfile,"(a7,I0,a4)") 'qd_hex_',int(mmax),".xyz"
	
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
	
	call atomcount(a,mmax,natom)
	
	write(300,*) natom
	write(300,*) "qdot"

	!call circle(rcut)
	!construindo a rede m

	!rede m parte de cima

	do n=0,mmax/2

		do m=0,mmax-1-n

			r(1)= dble(m)*r1(1)+dble(n)*r2(1)
			r(2)=dble(m)*r1(2)+dble(n)*r2(2)




				write(300,*) trim(tm),"      ",r(1),r(2),0.0

			

		end do


	end do

	!rede m parte de baixo

	do n=-1,-mmax/2+1,-1

		do m=abs(n),mmax-1

			r(1)= dble(m)*r1(1)+dble(n)*r2(1)
			r(2)=dble(m)*r1(2)+dble(n)*r2(2)


	

				write(300,*) trim(tm),"      ",r(1),r(2),0.0

			

		end do


	end do

	!construindo a rede x

	!rede x parte de cima

	do n=0,mmax/2-1

		do m=0,mmax-1-n

			r(1)= dble(m)*r1(1)+dble(n)*r2(1)+rp(1)
			r(2)=dble(m)*r1(2)+dble(n)*r2(2)+rp(2)


	


				write(300,*) trim(ch),"      ",r(1),r(2), (a/(sqrt(3.)*cos(theta)))*sin(theta)
				write(300,*) trim(ch),"      ",r(1),r(2),-(a/(sqrt(3.)*cos(theta)))*sin(theta)

			

		end do


	end do

	!rede x parte de baixo

	do n=-1,-mmax/2,-1

		do m=abs(n),mmax-1

			r(1)= dble(m)*r1(1)+dble(n)*r2(1)+rp(1)
			r(2)=dble(m)*r1(2)+dble(n)*r2(2)+rp(2)


	

				write(300,*) trim(ch),"      ",r(1),r(2), (a/(sqrt(3.)*cos(theta)))*sin(theta)
				write(300,*) trim(ch),"      ",r(1),r(2),-(a/(sqrt(3.)*cos(theta)))*sin(theta)

			

		end do


	end do


	!do n=-1,mmax

		!do m=0,mmax-n

			!r(1)= dble(m)*r1(1)+dble(n)*r2(1)+rp(1)
			!r(2)=dble(m)*r1(2)+dble(n)*r2(2)+rp(2)


	
	
				!write(201,*) m,n

				!write(300,*) "S","      ",r(1),r(2), (a/(sqrt(3.)*cos(theta)))*sin(theta)
				!write(300,*) "S","      ",r(1),r(2),-(a/(sqrt(3.)*cos(theta)))*sin(theta)


			

		!end do


	!end do


	close(300)


end program main

subroutine atomcount(a,mmax,natom)

	implicit none
	
	double precision :: a
	double precision,parameter :: theta=0.710

	double precision,dimension(2) :: r1,r2,rp,r
	integer :: m,n,i,j,erro,natom

	integer :: mmax !tamanho do ponto quantico
	double precision :: modv	

	r1(1)= a
	r1(2)= 0.

	r2(1)= a/2.
	r2(2)= (sqrt(3.)/2.)*a

	rp(1)= 0.
	rp(2)= a/sqrt(3.)
	
	natom = 0
	

	do n=0,mmax/2

		do m=0,mmax-1-n

			r(1)= dble(m)*r1(1)+dble(n)*r2(1)
			r(2)=dble(m)*r1(2)+dble(n)*r2(2)




				natom = natom + 1

			

		end do


	end do

	!rede m parte de baixo

	do n=-1,-mmax/2+1,-1

		do m=abs(n),mmax-1

			r(1)= dble(m)*r1(1)+dble(n)*r2(1)
			r(2)=dble(m)*r1(2)+dble(n)*r2(2)


	

				natom = natom + 1

			

		end do


	end do

	!construindo a rede x

	!rede x parte de cima

	do n=0,mmax/2-1

		do m=0,mmax-1-n

			r(1)= dble(m)*r1(1)+dble(n)*r2(1)+rp(1)
			r(2)=dble(m)*r1(2)+dble(n)*r2(2)+rp(2)


	


				natom = natom + 2

			

		end do


	end do

	!rede x parte de baixo

	do n=-1,-mmax/2,-1

		do m=abs(n),mmax-1

			r(1)= dble(m)*r1(1)+dble(n)*r2(1)+rp(1)
			r(2)=dble(m)*r1(2)+dble(n)*r2(2)+rp(2)


	

				natom = natom + 2

			

		end do


	end do	


end subroutine atomcount

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
