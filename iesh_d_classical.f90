program IESH_D_classical_2_dof

implicit none

!constant parameters
real*8, parameter :: clight=2.99792458D10,av=6.0221367D23,hbar=1.05457266D-34
real*8, parameter :: kb=1.3806503d-23
real*8, parameter :: pascal_to_atm=9.86923267d-6,kcal_to_J=6.95d-21
real*8, parameter :: amu2kg=1.66053892d-27
real*8, parameter ::au2kg=9.10938291d-31,au2J=4.35974417d-18,au2m=5.2917721092d-11,au2s=2.418884326505d-17
real*8, parameter :: q_el=2.307075d-28
complex*16,parameter :: iota = (0.d0,1.d0)
real*8 pi,wave_to_J,conv

!system variables
integer::n_el,nbasis,nclass,nmetal,ielec_hop
real*8,allocatable::H_diab(:,:),delH_dels(:,:,:),e_metal(:),mass(:),omg(:),Vc(:),H_diab_old(:,:)
real*8::band_width,gama_coup,omg1,omg2,g_coup1,g_coup2,V_reorg1,gamma_B,V_reorg2,omg1_prime,V_exothermicity
integer,allocatable::state(:),state_tentative(:)

!dynamical variables
real*8,allocatable::x(:),v(:),acc(:),acc_cl(:),acc_qm(:)
real*8,allocatable::x_old(:),v_old(:),acc_old(:)
real*8::dt,pot_cl,pot_en,pot_qm,t_max,ke_en
integer::nsteps,num_traj,traj_no

!other variables
real*8::temperature,beta,curr_time
integer::step_no,vib_state
real*8::t1,t2,classical_pot

!output variables
real*8,allocatable::pop(:,:),pop_gr(:,:),pop_ex(:,:),V_k(:)
integer::nstep_avg,time_steps,output_write,size_seed,seed2(2),i,seed_out
real*8::write_time,rnd,rnd_num
integer,allocatable::seed(:)

integer::ii
real*8::en,barrier

pi=4.d0*datan(1.d0)
wave_to_J=2*pi*clight*hbar
conv=1.d3/av

n_el=15
nclass=2
nmetal=30
nbasis= nmetal+1              !20 metal,2 molecular state
dt=0.25d-15
temperature=700.d0
beta=1.d0/(kb*temperature)
t_max=1000.d-12
nsteps=t_max/dt+1
num_traj=100
nstep_avg=100
time_steps=nsteps/nstep_avg
write_time=nstep_avg*dt
vib_state=12

allocate(H_diab(nbasis,nbasis),delH_dels(nbasis,nbasis,nclass),e_metal(nmetal),Vc(nmetal),H_diab_old(nbasis,nbasis))
allocate(mass(nclass),omg(nclass),state(n_el),state_tentative(n_el))
allocate(x(nclass),v(nclass),acc(nclass),x_old(nclass),v_old(nclass),acc_old(nclass),acc_cl(nclass),acc_qm(nclass))
allocate(pop(2,time_steps),pop_gr(vib_state,time_steps),pop_ex(vib_state,time_steps),V_k(nbasis))

call cpu_time(t1)
call system_parameters

open(10,file='ifolder.inp')
read(10,*)seed_out
close(10)
seed2(1)=1245
seed2(2)=3157

  call random_seed(size=size_seed)
  allocate(seed(size_seed))
  do i=1,size_seed/2
    seed(i)=seed2(1)*(2*i+i*i-7)*seed_out
  enddo
  do i=size_seed/2+1,size_seed
    seed(i)=seed2(2)*(i/2+34-i**3)*seed_out
  enddo
  call random_seed(put=seed)

pop=0.d0

do traj_no=1,num_traj
call initialize
output_write=1
call pop_calculation(1)

        do step_no=1,nsteps
                curr_time=step_no*dt
!ii=findloc(state,1,1)
!if(ii==0)then
!        en=0.5d0*mass(1)*omg1**2*(x(1)-g_coup1/(mass(1)*omg1**2))**2+0.5d0*mass(1)*omg2**2*(x(2)-g_coup2*x(1)/(mass(1)*omg2**2))**2
!else
!        en=0.5d0*mass(1)*omg1**2*(x(1)+g_coup1/(mass(1)*omg1**2))**2+0.5d0*mass(1)*omg2**2*(x(2)-g_coup2*x(1)/(mass(1)*omg2**2))**2+V_exothermicity
!endif
!if(mod(step_no,1000)==0) write(300,*)curr_time*1.d12,x*1.d10

!write(201,'(1es14.7,20i)')curr_time*1.d12,state
!write(202,'(25es14.7)')curr_time*1.d12,V_k(state)/wave_to_J

if(step_no>00)                call LZ_hop
!                if(step_no==1) call pop_calculation(1)
                if(mod(step_no,nstep_avg)==0)then
                output_write=output_write+1
                call pop_calculation(output_write)
                endif
                call evolution
!                call output
        enddo
enddo

call write_output

call cpu_time(t2)
!write(15,*)t2-t1

contains
!---------------------------------------------------------------------------------------------------------------------------
subroutine system_parameters
implicit none
integer::i,j
real*8::knots(nmetal/2),weights(nmetal/2),VER_rate,w

mass=1836.d0*au2kg
omg(1)=2000.d0*(2*pi*clight)
omg(2)=400.d0*(2*pi*clight)
omg1=omg(1)
omg2=omg(2)
VER_rate=0.1d12
V_reorg1=wave_to_J*20000.d0
gamma_B=2*pi*clight*400.d0
band_width=0.064*au2J
gama_coup=au2J*1.d-4!*4.d0
V_exothermicity=wave_to_J*(5000.d0)

w=omg1
g_coup1=dsqrt(V_reorg1*mass(1)*omg1**2/2.d0)
V_reorg2=VER_rate*mass(1)*omg1
V_reorg2=V_reorg2 * 2 *((w*w-omg2*omg2)**2+(gamma_B*w)**2)/(omg2**2*gamma_B*w)
g_coup2=dsqrt(V_reorg2*mass(1)*omg2**2/2.d0)

omg1_prime=sqrt(omg1**2+g_coup2**2/(mass(1)**2*omg2**2))

!write(6,*)-g_coup1/(mass(1)*omg1**2)*1.d10
!write(6,*)(-2*g_coup1**2/(mass(1)*omg1**2)+V_exothermicity)/wave_to_J
!stop


open(30,file="15_knot_points_x.txt")
!open(30,file="knot_points_x.inp")
  do i=1,nmetal/2
    read(30,*) knots(i)
  enddo
  close(30)

open(30,file="15_knot_points_w.txt")
!open(30,file="knot_points_w.inp")
  do i=1,nmetal/2
    read(30,*) weights(i)
  enddo
  close(30)

  do i=1,nmetal/2
    e_metal(nmetal/2-i+1)=-band_width/2.d0*(1/2.d0+knots(i)/2.d0)
    e_metal(nmetal/2+i)=band_width/2.d0*(1/2.d0+knots(i)/2.d0)
    Vc(nmetal/2-i+1)=sqrt(gama_coup/(2*pi))*sqrt(band_width*weights(i))/2.d0
    Vc(nmetal/2+i)=sqrt(gama_coup/(2*pi))*sqrt(band_width*weights(i))/2.d0
  enddo


end subroutine
!-----------------------------------------------------------------------------------------------------------------------------------------------------------
subroutine initialize
implicit none
integer::i,init_vib_state,j
real*8::sig_x,sig_p,rnd,e_vib,theta,en,omg2_prime,x2_prime,barrier

init_vib_state=10
state=0
!state(1)=1
do i=1,n_el
state(i)=i+1
enddo

do i=1,nclass
    sig_x=1.d0/dsqrt(beta*mass(i)*omg(i)**2)
    sig_p=dsqrt(mass(i)/beta)
    call gaussian_random_number(rnd)
    x(i)=rnd*sig_x
    call gaussian_random_number(rnd)
    v(i)=(1.d0/mass(i)*(rnd*sig_p))
enddo

x(1)=x(1)-g_coup1/(mass(1)*omg1**2)
x(2)=x(2)-g_coup2*x(1)/(mass(1)*omg2**2)

e_vib=(init_vib_state-0.5d0)*hbar*omg1
call random_number(rnd)
theta=rnd*2*pi
!omg2_prime=dsqrt((omg(2)**2-g_coup2**2/(mass(1)**2*omg1_prime**2)))
!x2_prime=0.5d0*g_coup1*g_coup2/(mass(1)**2*omg1_prime**2*omg2_prime**2)

!x(2)=x(2)+x2_prime!+g_coup2*x(1)/(mass(1)*omg2**2)
!x(1)=dsqrt(2*e_vib/(mass(1)*omg1_prime**2))*dcos(theta)+(g_coup1+g_coup2*x(2))/(mass(1)*omg1_prime**2)
x(1)=dsqrt(2*e_vib/(mass(1)*omg1**2))*dcos(theta)+g_coup1/(mass(1)*omg1**2)
v(1)=dsqrt(2*e_vib/mass(1))*dsin(theta)
x(2)=x(2)-g_coup2*x(1)/(mass(1)*omg2**2)


!en=0.5d0*mass(1)*omg1**2*(x(1)-g_coup1/(mass(1)*omg1**2))**2+0.5d0*mass(1)*v(1)**2
!en=en/(hbar*omg1)
!write(6,*)en,x*1.d10
!enddo
!stop

call TDSE

end subroutine
!-----------------------------------------------------------------------------------------------------------------------------------------------------------
subroutine evolution
implicit none
real*8::gama_dt,c1,c2,c0,delta_r(nclass),delta_v(nclass)

H_diab_old=H_diab

    x(1)=x(1)+v(1)*dt+0.5*acc(1)*dt*dt
!    x(2)=x(2)+v(2)*dt+0.5*acc(2)*dt*dt
    gama_dt=gamma_B*dt
    c0=dexp(-gama_dt)
    c1=1.d0/gama_dt*(1.d0-c0)
    c2=1.d0/gama_dt*(1.d0-c1)
    call stochastic_force(delta_r,delta_v,dt)
    x(2)=x(2)+c1*dt*v(2)+c2*dt*dt*acc(2)+delta_r(2)
    v(1)=v(1)+0.5*acc(1)*dt
!    v(2)=v(2)+0.5*acc(2)*dt
    acc_old=acc
    call TDSE
    v(1)=v(1)+0.5*dt*acc(1)
!    v(2)=v(2)+0.5*acc(2)*dt
    v(2)=c0*v(2)+(c1-c2)*dt*acc_old(2)+c2*dt*acc(2)+delta_v(2)

!call TDSE

end subroutine
!-----------------------------------------------------------------------------------------------------------------------------------------------------------
subroutine equilibration
implicit none

integer::i,j,k

call eq_acc

do i=1,10000
x=x+v*dt+0.5d0*acc_cl*dt*dt
v=v+0.5d0*acc_cl*dt
call eq_acc
v=v+0.5d0*acc_cl*dt
write(100,*)i,x*1.d10
write(101,*)i,v
enddo


end subroutine 
!-----------------------------------------------------------------------------------------------------------------------------------------------------------
subroutine LZ_hop
implicit none
integer::i,j,k,ii
real*8::p_LZ,rnd,delF,V1,V2,en_gr,en_ex

ii=0
ii=findloc(state,1,1)

if(ii==0)then      !neutral state active
       outer: do i=1,n_el
                V1=H_diab(1,1)-H_diab(state(i),state(i))!-H_diab(2,2)
                V2=H_diab_old(1,1)-H_diab_old(state(i),state(i))!-H_diab_old(2,2)
                if(V1*V2<0.d0)then
                        delF=2.d0*g_coup1
                        p_LZ=abs(2*pi*Vc(state(i)-1)**2/(hbar*v(1)*delF))
                call random_number(rnd)
                if(rnd<p_LZ) then
!        en_ex=0.5d0*mass(1)*omg1**2*(x(1)+g_coup1/(mass(1)*omg1**2))**2+0.5d0*mass(1)*v(1)**2
!        en_ex=en_ex/(hbar*omg1)
!        en_gr=0.5d0*mass(1)*omg1**2*(x(1)-g_coup1/(mass(1)*omg1**2))**2+0.5d0*mass(1)*v(1)**2
!        en_gr=en_gr/(hbar*omg1)
!        write(200,'(2i,5es14.7)')floor(en_gr),floor(en_ex),0.d0,(H_diab(1,1))/wave_to_J,H_diab(state(i),state(i))/wave_to_J
                        state(i)=1
                        exit outer
                endif
                endif
        enddo outer
  else
outer1: do j=2,nbasis
      if(.not.(any(j==state))) then
          V1=H_diab(1,1)-H_diab(j,j)!-H_diab(2,2)
          V2=H_diab_old(1,1)-H_diab_old(j,j)!-H_diab_old(2,2)
          if(V1*V2<0.d0) then
                delF = 2.d0*g_coup1!2*g_coup2 * g_coup1/(mass(1)*omg1_prime**2)
                p_LZ=abs(2*pi*Vc(j-1)**2/(hbar*v(1)*abs(delF)))! * fc**2
                call random_number(rnd)
                if(rnd<p_LZ) then
!        en_ex=0.5d0*mass(1)*omg1**2*(x(1)+g_coup1/(mass(1)*omg1**2))**2+0.5d0*mass(1)*v(1)**2
!        en_ex=en_ex/(hbar*omg1)
!        en_gr=0.5d0*mass(1)*omg1**2*(x(1)-g_coup1/(mass(1)*omg1**2))**2+0.5d0*mass(1)*v(1)**2
!        en_gr=en_gr/(hbar*omg1)
!        write(200,'(2i,5es14.7)')floor(en_gr),floor(en_ex),1.d0,(H_diab(1,1))/wave_to_J,H_diab(j,j)/wave_to_J
        !write(200,*)'a to n',floor(en_gr),floor(en_ex)
                        state(ii)=j
                        exit outer1
                endif
          endif
      endif
      enddo outer1
endif
!call TDSE
end subroutine

!-----------------------------------------------------------------------------------------------------------------------------------------------------------
subroutine pop_calculation(j)
implicit none
integer,intent(in)::j
integer::i,ii
real*8::en

ii=0
ii=findloc(state,1,1)

en=0.d0
if(ii.ne.0)then
        pop(2,j)=pop(2,j)+1.d0
        !en=0.5d0*mass(1)*omg1_prime**2*(x(1)+(g_coup1-g_coup2*x(2))/(mass(1)*omg1_prime**2))**2+0.5d0*mass(1)*v(1)**2
        en=0.5d0*mass(1)*omg1**2*(x(1)+g_coup1/(mass(1)*omg1**2))**2+0.5d0*mass(1)*v(1)**2
        en=en/(hbar*omg1)
        do i=1,vib_state
        if(en>(i-1.d0).and.en<real(i)) pop_ex(i,j)=pop_ex(i,j)+1.d0
        enddo
else
        pop(1,j)=pop(1,j)+1.d0
        !en=0.5d0*mass(1)*omg1_prime**2*(x(1)-(g_coup1+g_coup2*x(2))/(mass(1)*omg1_prime**2))**2+0.5d0*mass(1)*v(1)**2
        en=0.5d0*mass(1)*omg1**2*(x(1)-g_coup1/(mass(1)*omg1**2))**2+0.5d0*mass(1)*v(1)**2
        en=en/(hbar*omg1)
        do i=1,vib_state
        if(en>(i-1.d0).and.en<real(i)) pop_gr(i,j)=pop_gr(i,j)+1.d0
        enddo
endif


do i=1,nbasis
V_k(i)=H_diab(i,i)
enddo

!do i=1,n_el
!if(state(i)==1)then
!pop(2,j)=pop(2,j)+1.d0
!endif
!enddo

end subroutine
!-----------------------------------------------------------------------------------------------------------------------------------------------------------
subroutine write_output
implicit none
integer::i

pop=pop/real(num_traj)
pop_gr=pop_gr/real(num_traj)
pop_ex=pop_ex/real(num_traj)

open(10,file='pop.out')
open(11,file='pop_vib.out')
!open(12,file='pop_ex_1.out')

do i=1,time_steps
write(10,*)i*write_time*1.d12,pop(1,i),pop(2,i)
write(11,'(40es14.5)')i*write_time*1.d12,pop_gr(:,i),pop_ex(:,i)
!write(12,'(20es14.5)')i*write_time*1.d12,pop_ex(:,i)
enddo

end subroutine
!-----------------------------------------------------------------------------------------------------------------------------------------------------------
subroutine TDSE
implicit none
integer::i,j,k
real*8::pot1,pot2

H_diab=0.d0
delH_dels=0.d0

!H_diab(1,1)=0.5d0*mass(1)*omg1_prime**2*(x(1)-(g_coup1-g_coup2*x(2))/(mass(1)*omg1_prime**2))**2-0.5d0*(g_coup1-g_coup2*x(2))**2/(mass(1)*omg1_prime**2)
!H_diab(2,2)=0.5d0*mass(1)*omg1_prime**2*(x(1)+(g_coup1+g_coup2*x(2))/(mass(1)*omg1_prime**2))**2-0.5d0*(g_coup1+g_coup2*x(2))**2/(mass(1)*omg1_prime**2)+V_exothermicity

!H_diab(1,1)=0.5d0*mass(1)*omg(1)**2*(x(1)-g_coup1/(mass(1)*omg(1)**2))**2
!H_diab(2,2)=0.5d0*mass(1)*omg(1)**2*(x(1)+g_coup1/(mass(1)*omg(1)**2))**2+V_exothermicity

H_diab(1,1)=2*g_coup1*x(1)+V_exothermicity
delH_dels(1,1,1)=2*g_coup1

do i=2,nbasis
H_diab(i,i)=e_metal(i-1)
enddo

!delH_dels(1,1,1)=mass(1)*omg1_prime**2*(x(1)-(g_coup1-g_coup2*x(2))/(mass(1)*omg1_prime**2))
!delH_dels(1,1,2)=g_coup2*x(1)
!delH_dels(2,2,1)=mass(1)*omg1_prime**2*(x(1)+(g_coup1+g_coup2*x(2))/(mass(1)*omg1_prime**2))
!delH_dels(2,2,2)=g_coup2*x(1)

!delH_dels(1,1,1)=mass(1)*omg(1)**2*x(1)-g_coup1
!delH_dels(1,1,2)=0.d0
!delH_dels(2,2,1)=mass(1)*omg(1)**2*x(1)+g_coup1
!delH_Dels(2,2,2)=0.d0

acc_cl=0.d0
!acc_cl(1)=-g_coup2*(x(2)-g_coup2*x(1)/(mass(1)*omg(2)**2))
!acc_cl(2)=mass(1)*omg(2)**2*(x(2)-g_coup2*x(1)/(mass(1)*omg(2)**2))
!acc_cl=-acc_cl/mass

acc_cl(1)=mass(1)*omg(1)**2*x(1)-g_coup1-g_coup2*x(2)+g_coup2**2*x(1)/(mass(1)*omg(2)**2)
acc_cl(2)=mass(1)*omg(2)**2*(x(2)-g_coup2*x(1)/(mass(1)*omg(2)**2))

acc_cl=-acc_cl/mass
acc_qm=0.d0
do i=1,nclass
do j=1,n_el
if(state(j)>0) acc_qm(i)=acc_qm(i)-delH_dels(state(j),state(j),i)/mass(i)
enddo
enddo

!do i=1,nmetal
!write(6,*)H_diab(i+2,i+2)/conv
!enddo
!write(6,*)H_diab(1,1)/conv,H_diab(2,2)/conv
!stop

acc=acc_cl+acc_qm
end subroutine
!---------------------------------------------------------------------------------------------------------------------------------------------------------------
subroutine eq_acc
implicit none

integer::i,j

acc_cl=0.d0
acc_cl(1)=mass(1)*omg(1)**2*x(1)-g_coup1-g_coup2*x(2)+g_coup2**2*x(1)/(mass(1)*omg(2)**2)
acc_cl(2)=mass(1)*omg(2)**2*(x(2)-g_coup2*x(1)/(mass(1)*omg(2)**2))

acc_cl=-acc_cl/mass
end subroutine
!---------------------------------------------------------------------------------------------------------------------------------------------------------------
subroutine output
implicit none
integer::i

pot_en=0.d0
do i=1,n_el
pot_en=pot_en+H_diab(state(i),state(i))
enddo
pot_en=pot_en+0.5d0*mass(1)*omg(1)**2*(x(1)-g_coup1/(mass(1)*omg(1)**2))**2+0.5d0*mass(1)*omg(2)**2*(x(2)-g_coup2*x(1)/(mass(1)*omg(2)**2))**2

ke_en=0.d0
do i=1,nclass
ke_en=ke_en+0.5d0*mass(i)*v(i)**2
enddo

write(50,*)curr_time*1.d12,pot_en/wave_to_J,ke_en/wave_to_J
write(51,*)curr_time*1.d12,x*1.d10
end subroutine
!---------------------------------------------------------------------------------------------------------------------------------------------------------------
subroutine stochastic_force(delr,delv,dt)
  implicit none
  real*8,intent(in)::dt
  real*8,intent(out) :: delr(nclass),delv(nclass)!f(nclass)
  integer i
  real*8 rnd1,rnd2,sig_r,sig_v,sig_rv,gdt

  gdt=gamma_B*dt
  do i=1,nclass
    sig_r=dt*dsqrt(kb*temperature/mass(i)*1.d0/gdt*(2-1.d0/gdt*(3-4*dexp(-gdt)+dexp(-2*gdt))))
    sig_v=dsqrt(kb*temperature/mass(i)*(1-dexp(-2*gdt)))
    sig_rv=(dt*kb*temperature/mass(i)* 1.d0/gdt*(1-dexp(-gdt))**2)/(sig_r*sig_v)  !! correlation coeffecient

    call gaussian_random_number(rnd1)
    call gaussian_random_number(rnd2)
    delr(i)=sig_r*rnd1
    delv(i)=sig_v*(sig_rv*rnd1+dsqrt(1-sig_rv**2)*rnd2)
  enddo

!  delr=delr-sum(delr)/dfloat(nclass)
!  delv=delv-sum(delv)/dfloat(nclass)

end subroutine stochastic_force
!-----------------------------------------------------------------
subroutine gaussian_random_number(rnd)
  !! generates gaussian distribution with center 0, sigma 1
  !! q0+sig*rnd gives center=q0, sigma=sig
  implicit none
  real*8,intent(out)::rnd
  real*8 rnd1,rnd2,pi

  pi=dacos(-1.d0)
  call random_number(rnd1)
  call random_number(rnd2)
  rnd = dsqrt(-2*log(rnd1))*dcos(2*pi*rnd2)

end subroutine gaussian_random_number
!-----------------------------------------------------------------
end program

