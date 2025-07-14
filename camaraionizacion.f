!Programa de simulación Monte Carlo de una cámara de ionización esférica 
program camaraionizacion
use randomnumber
implicit none
real*8::h,PI,H_elec,rhomed,rho,Vol,hmax,hmin,haux,dt,alfa,mu_pos,mu_neg,k,t,V,recol_neg,recol_pos,a,b,f,fsum,fcuadsum,error,fmed,er
real*8,dimension(:),allocatable::rholocal,r_neg,r_pos,phi_neg,phi_pos,theta_neg,theta_pos,p
integer::N0,Npos,Nneg,i,j,m,i3,i4,m1

!se calcula el volumen interior de la esfera conociendo los radios de los electrodos exterior(a) e interior(b)
a=0.305
b=0.055
PI=3.141592
Vol=4/3.0*PI*(a**3-b**3)

!Calculamos la densidad de iones sabiendo que una radiación de 1 mGy genera 2.183E8/cm³ y dividirla por un cierto factor que luego se añadirá a la constante de recombinación
rho=2.183E8/(5.0E05)
!Tomamos como constante de recombinación k =1.6E-6 cm³/s  y multiplicamos por el mismo factor que habiamos dividido la densidad
k=(1.6E-6)*(5.0E5)
!Calculamos el número inicial de iones
N0=int(rho*Vol)
Nneg=N0
Npos=NNeg

!Debemos elegir un voltaje, así como un valor a las movilidades en cm²/(V*s) y el parámetro alfa.
V=50
mu_pos=1.80
mu_neg=2.10
alfa=100

allocate(rholocal(1:Nneg))
allocate(r_neg(1:Nneg))
allocate(r_pos(1:Npos))
allocate(phi_neg(1:Nneg))
allocate(phi_pos(1:Npos))
allocate(theta_neg(1:Nneg))
allocate(theta_pos(1:Npos))
allocate(p(1:Nneg))

!inicializando algunas variables
error=1
er=1
m1=1
fsum=0.0
fcuadsum=0.0
!Bucle buscando error relativo pequeño
do while(er.gt.(0.001))
	
	!inicializamos de nuevo el numero de iones(que antes ha disminuido), así como el tiempo y la recolección.
	Nneg=N0
	Npos=NNeg
	t=0.0
	recol_pos=0
	recol_neg=0
	
	!Podríamos crear un GIF para un corte de la esfera de un cierto grosor
	!OPEN(unit=1,file='posicionesneg',status='old')
	!OPEN(unit=2,file='posicionespos',status='old')
	
	!Generamos una distribucion inicial de iones.
	call dran_ini(m1*1000+111813)	
	call distribucionini
	!hallamos un valor apropiado de h.
	call calcularh
	!Calculamos el intervalo temporal conociendo el campo eléctrico máximo, el cual será para el radio mínimo->b
	dt=h/(mu_neg*E(b)*alfa)
	!Bucle de pasos temporales de cada simulación que terminará cuando no haya iones que pueda recombinarse
	do while (rhomed.gt.0)
		!Vamos a escribir la posicion de cada ion negativo y positivo con objetivo de hacer gift
!		do i=1,N0
!			if((i.le.Nneg).and.(r_neg(i)*cos(theta_neg(i)).ge.(-0.05)).and.(r_neg(i)*cos(theta_neg(i)).le.(0.05)))then
			!write(1,*)r_neg(i)*cos(phi_neg(i))*sin(theta_neg(i)),r_neg(i)*sin(phi_neg(i))*sin(theta_neg(i)),r_neg(i)*cos(theta_neg(i))
!			else 
!			write(1,*) 10,10,10
!			end if
!			end do
!			write(1,*)' '
!			write(1,*)' '
!			do i=1,N0		
!			if((i.le.Npos).and.(r_pos(i)*cos(theta_pos(i)).ge.(-0.05)).and.(r_pos(i)*cos(theta_pos(i)).le.(0.05)))then
			!write(2,*)r_pos(i)*cos(phi_pos(i))*sin(theta_pos(i)),r_pos(i)*sin(phi_pos(i))*sin(theta_pos(i)),r_pos(i)*cos(theta_pos(i))
!			else 
!			write(2,*) 10,10,10
!			end if
!			end do
!			write(2,*)' '
!			write(2,*)' '
		!Hacemos evolucionar el sistema
		call pasotemporal
		
	end do
	!Cuando finaliza una simulación calculamos la eficiencia y su incertidumbre
	f=100.0*(recol_neg+Nneg)/(1.0*N0)
	fsum=fsum+f
	fcuadsum=fcuadsum+f**2
	fmed=fsum/(1.0*m1)
	!Obligar al programa a hacer un mínimo de 10 simulaciones
	if(m1.gt.9)then
		error=sqrt(1/(m1*1.0*(m1-1))*(fcuadsum-fsum**2/(1.0*m1)))
		er=error/fmed
	end if
	!Enviar resultados a la pantalla
	write(6,*)m1,' ',fmed,' ',error,' ',er
	m1=m1+1

end do

stop
contains
!Aquí encontramos las distintas subrutinas y funciones empleadas.

!Función para calcular la distancia entre un ión negativo 0 y un ion positivo l
real*8 function dist(o,l)
	integer::o,l
	real*8::auxx,auxy,auxz
	auxx=r_pos(l)*cos(phi_pos(l))*sin(theta_pos(l))-r_neg(o)*cos(phi_neg(o))*sin(theta_neg(o))
	auxy=r_pos(l)*sin(phi_pos(l))*sin(theta_pos(l))-r_neg(o)*sin(phi_neg(o))*sin(theta_neg(o))
	auxz=r_pos(l)*cos(theta_pos(l))-r_neg(o)*cos(theta_neg(o))
	dist=sqrt(auxx**2+auxy**2+auxz**2)
end function
		
!Función para calcular el módulo del campo eléctrico en un punto en funcion de diferencia de potencial entre electrodos
real*8 function E(r1)
	real*8::r1
	E=V*a*b/((r1**2)*(a-b))
end function

!Subrutina para hallar h mediante el algoritmo de bisección	
subroutine calcularh
	real*8::aux
	!Elegimos un valor inicial del error(aux) mayor que 0.01 para que se produzca el do while
	aux=1
	!Implementamos el algoritmo de biseccion estimando un intervalo que puede contener a h, por ejemplo el semiradio mayor
	hmax=a
	hmin=0.0
	!Realizamos el bucle las veces que sea necesario hasta alcanzar un valor apropiado de h. Esto se hará cuando el error relativo sea menor del 1%.
	do while(aux.gt.(0.01))
		!do m=1,3
		h=(hmax+hmin)/2.0
		!write(6,*) 'El valor de h es: ',h
		call densidades
		!veamos el error relativo de pholocal para un h dado
		aux=abs(rhomed-rho)/rho
		!write(6,*)'Densidad promedio de densidades locales: ',rhomed
		!write(6,*)'Densidad total: ',rho
		!write(6,*)'error relativo', aux
		if(rhomed.gt.rho) then
			!hmax=h
			hmin=h
		else 
			!hmin=h
			hmax=h
		end if
	end do
	!end do		
end subroutine

!Esta subrutina que crea una distribución uniforme a partir de un cierto número de iones.
subroutine distribucionini
	do m=1,Nneg
		!El volumen de una esfera crece con r³ luego para subsanar este efecto le calcularemos la raiz cúbica a la distribucion unitaria aleatoria haciendo que ahora sí la densidad se mantenga constante en todo el disco y no haya acumulación de partículas para radios pequeños.
		r_neg(m)=b+(a-b)*(abs(dran_u()))**(1.0/3)
		r_pos=r_neg
		phi_neg(m)=abs(dran_u()*2*PI)
		phi_pos(m)=phi_neg(m)
		theta_neg(m)=abs(dran_u()*PI)
		theta_pos(m)=theta_neg(m)
	end do
end subroutine

!Esta subrutina permite calcular las densidades(positivas) locales asociadas a cada ion negativo, así como la densidad local media. 
subroutine densidades
	rholocal=0
	rhomed=0
		!Primero debemos hallar la densidad local positiva asociada a cada ion  negativo para poder realizar el promedio, para ello sumaremos la contribucion de los iones positivos j que rodean a un ion negativo i a una distancia menor que h.
		do i=1,Nneg
			do j=1,Npos
				if(dist(i,j).lt.h)then
					rholocal(i)=rholocal(i)+3.0/(4*PI*h**3) 
					!write(6,*) rholocal(i)
					!write(6,*)'Densidad local de ion numero ',i,' :',rholocal(i)
					!write(6,*)'El ionnegativo nº ',i,'interactua con el ion positivo ',j,' por proximidad'
				end if
			
			end do	
			rhomed=rhomed+rholocal(i)/(Nneg*1.0)
		end do
end subroutine


!Esta subrutina permite evolucionar un sistema dado desde un instante t a un instante t+dt
subroutine pasotemporal
	integer::l
	real*8::aux
	!En primer lugar elegimos el valor del salto temporal dt.Para ello emplearemos la mobilidad iónica negtiva y el h hallado anteriormente. El dt será el tiempo que tarda una partícula en moverse 				una distancia h, estando sometida la particula a un campo electrico E.
	
	call densidades
	!En segundo lugar tenemos en cuenta la recombinación.
	l=1
	do while (l.le.Nneg)
		!Calculamos la probabilidad de recombinación aproximada o exacta
		!p(l)=k*rholocal(l)*dt
		p(l)=1-exp(-k*rholocal(l)*dt)
		!write(6,*)'La probabilidad de recombinación del ion negativo ',l,' es',p(l)
		aux=dran_u()
		!Eliminamos el ion negativo si la probabilidad es mayor que el nº aleatorio obtenido. Tmb eliminamos el ion positivo mas cercano
		if(p(l).gt.aux)then
			call eliminariones(l)
			!reducimos el indice l para que se repita el calculo de probabilidad para la partícula N, ahora en el lugar l al haber elminiado la partícula l 
			l=l-1
		end if
		l=l+1
	end do
	
	!Por ultimo movemos cada partícula, primero las negativas y luego las positivas.
	do l=1,Nneg
		r_neg(l)=r_neg(l)-mu_neg*E(r_neg(l))*dt
		
		!veamos que iones son recolectados
		!if(y_neg(l).lt.0)then
		!	recol_neg=recol_neg+1
		!	Write(6,*)'El ion negativo ',l,' ha alcanzado el electrodo positivo'
			
		!end if
	end do
	do l=1,Npos
		r_pos(l)=r_pos(l)+mu_pos*E(r_pos(l))*dt
		!if (y_pos(l).gt.d)then
		!	recol_pos=recol_pos+1
		!	Write(6,*)'El ion positivo ',l,' ha alcanzado el electrodo negativo'
		!end if
	end do
	call reordenar
	t=t+dt
end subroutine

!Esta subrutina elimina el ión negativo i2 y el ión positivo más cercano.
subroutine eliminariones(i2)
integer::iaux,i2
iaux=1
	!Veamos cual es el ion positivo m mas cercano al negativo i2. Empezamos el bucle comparando las distancias a i de lios iones positivos 1 y 2
	do m=2,Npos
		if(dist(i2,m).lt.dist(i2,iaux))then
			iaux=m
			
		end if
	end do
	!Como realizamos la eliminación de iones? Por ejemplo sustituimos el N por el ion i o k en cuestion.
	!write(6,*)' Los iones negativo nº',i2,' y positivo nº ',iaux,' se han recombinado.'
	
	r_neg(i2)=r_neg(Nneg)
	phi_neg(i2)=phi_neg(Nneg)
	theta_neg(i2)=theta_neg(Nneg)
	rholocal(i2)=rholocal(Nneg)
	Nneg=Nneg-1	
	
	r_pos(iaux)=r_pos(Npos)
	phi_pos(iaux)=phi_pos(Npos)
	theta_pos(iaux)=theta_pos(Npos)
	Npos=Npos-1
		
end subroutine

!Esta subrutina reordena los coeficientes de los iones eliminando aquellos que se hayan recombinado o esten en los electrodos.
subroutine reordenar
	i4=1
	!ordenamos iones negativos
	do while(i4.le.Nneg)
		if(r_neg(i4).le.b)then
			r_neg(i4)=r_neg(Nneg)
			phi_neg(i4)=phi_neg(Nneg)
			theta_neg(i4)=theta_neg(Nneg)
			rholocal(i4)=rholocal(Nneg)
			Nneg=Nneg-1	
			recol_neg=recol_neg+1
			!write(6,*)'El ion negativo ',i4,' ha desaparecido.'	
		else
		i4=i4+1
		end if		
		
	end do
	i4=1	
	do while(i4.le.Npos)
		if(r_pos(i4).ge.a)then
			r_pos(i4)=r_pos(Npos)
			phi_pos(i4)=phi_pos(Npos)
			theta_pos(i4)=theta_pos(Npos)
			Npos=Npos-1
			recol_pos=recol_pos+1
			!write(6,*)'El ion positivo ',i4,' ha desaparecido.'		
		else
		i4=i4+1
		end if
		
	end do
end subroutine

end program
