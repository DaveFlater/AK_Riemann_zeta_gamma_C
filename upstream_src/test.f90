    program test
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    use gamma_zeta_module
    implicit none   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~        
    integer, parameter              :: qp=selected_real_kind(33, 4931)
    real(kind=qp), parameter        :: pi=4*atan(1.0_qp)
    complex(kind=qp), parameter     :: ii=(0.0_qp,1.0_qp)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    complex(kind=qp)        :: z, s, f(1:10), g(1:10), h(1:2)
    real(kind=qp)           :: error, t(1:5)
    real                    :: U(1:2)   
    integer                 :: k
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    print *, 'Test 1 : compute Gamma(k) for integer values of k'
    do k=1,10
        z=k
        f(k)=exp(ln_gamma(z))
        print *,f(k)
    end do
    g=(/1,1,2,6,24,120,720,5040, 40320,362880/)
    print *, 'max relative error=',maxval(abs(f/g-1))
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    print *, 'Test 2 : compute Gamma(1/2+k) for integer values of k'
    do k=1,5
        z=k-0.5_qp
        f(k)=exp(ln_gamma(z))
        print *,f(k)
    end do
    g(1:5)=sqrt(pi)*(/1.0_qp,0.5_qp,0.75_qp,1.875_qp,6.5625_qp/)
    print *, 'max relative error=',maxval(abs(f(1:5)/g(1:5)-1))
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    print *, 'Test 3 : compute |log(Gamma(z+1))-log(Gamma(z))-log(z)| for 10000 random values of z, which are'
    print *, '    uniformly distributed in the square |Im(z)|<20, |Re(z)|<20'
    error=0.0_qp
    do k=1,10000 
        call random_number(U)
        z=20*complex(2*U(1)-1,2*U(2)-1) 
        error=max(error,abs(ln_gamma(z+1)-ln_gamma(z)-log(z)))
    end do
    print *,'maximum error=',error  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    print *, 'Test 4 : check that psi(1)=-0.57721566490153286060651209008240243 (Euler-Mascheroni constant)'
    print *, 'psi(1)=',psi(1.0_qp+0*ii)
    print *, 'error=',psi(1.0_qp+0*ii)+0.57721566490153286060651209008240243_qp
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    print *, 'Test 5 : compute |psi(z+1)-psi(z)-1/z| for 10000 random values of z, which are'
    print *, '    uniformly distributed in the square |Im(z)|<20, |Re(z)|<20'
    error=0.0_qp
    do k=1,10000 
        call random_number(U)
        z=20*complex(2*U(1)-1,2*U(2)-1) 
        error=max(error,abs(psi(z+1)-psi(z)-1/z))
    end do
    print *,'maximum error=',error      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    print *, 'Test 6 : compute zeta(s) for s=2,4,6,8,10'
    error=0.0_qp
    do k=1,5 
        s=2*k
        f(k)=Riemann_zeta(s)
        print *,f(k)
    end do
    g(1)=pi**2/6
    g(2)=pi**4/90
    g(3)=pi**6/945
    g(4)=pi**8/9450
    g(5)=pi**10/93555    
    print *,'maximum error=',maxval(abs(f(1:5)-g(1:5))) 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    print *, 'Test 7 : compute zeta(s) for s=0,-1,-3,-5,-7,-9'
    error=0.0_qp
    s=0.0_qp    
    f(1)=Riemann_zeta(s)
    print *,f(1)
    do k=1,5 
        s=1-2*k
        f(k+1)=Riemann_zeta(s)
        print *,f(k+1)
    end do
    g(1)=-0.5_qp
    g(2)=-1.0_qp/12
    g(3)=1.0_qp/120
    g(4)=-1.0_qp/252
    g(5)=1.0_qp/240
    g(6)=-1.0_qp/132
    print *,'maximum error=',maxval(abs(f(1:6)-g(1:6)))
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! the zeros of zeta(s) used in the next three tests
! were downloaded from https://www-users.cse.umn.edu/~odlyzko/zeta_tables/index.html
    print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    print *, 'Test 8 : compute zeta(1/2+i*t_k) for the first five non-trivial zeros 1/2+i*t_k of zeta(s)'
    t(1)=14.13472514173469379045725198356247027_qp
    t(2)=21.02203963877155499262847959389690277_qp
    t(3)=25.01085758014568876321379099256282181_qp
    t(4)=30.42487612585951321031189753058409132_qp
    t(5)=32.93506158773918969066236896407490348_qp
    do k=1,5 
        print *, Riemann_zeta(0.5_qp+ii*t(k))
    end do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    print *, 'Test 9 : compute zeta(1/2+i*t_k) for the first five non-trivial zeros 1/2+i*t_k of zeta(s) at height t=200'
    t(1)=201.2647519437037887330161334275481732_qp
    t(2)=202.4935945141405342776866606378643158_qp
    t(3)=204.1896718031045543307164383863136851_qp
    t(4)=205.3946972021632860252123793906930909_qp
    t(5)=207.9062588878062098615019679077536442_qp
    do k=1,5 
        print *, Riemann_zeta(0.5_qp+ii*t(k))
    end do  
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    print *, 'Test 10 : compute zeta(1/2+i*t_k) for the first five non-trivial zeros 1/2+i*t_k of zeta(s) at height t=10^6'
    print *, '   note that these zeros are only correct to 10^{-9},'
    print *, '   so the values of zeta(1/2+i*t_k) should be of the same order of magnitude' 
    t(1)=1000000.584097696_qp
    t(2)=1000000.828343490_qp
    t(3)=1000001.435265267_qp
    t(4)=1000001.905648406_qp
    t(5)=1000002.877617740_qp
    do k=1,5 
        print *, Riemann_zeta(0.5_qp+ii*t(k))
    end do      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    print *, 'Test 11 : computing the first non-trivial zero with Newtons method'
    z=0.5_qp+14.13472514173469379045725198356247027_qp*ii 
    s=0.5_qp+14.1_qp*ii
    do k=1,6
        h=Riemann_zeta_prime(s)
        s=s-h(1)/h(2)
        print *,'next iteration of Newtons method:', s
    end do
    print *,'error:', abs(s-z)
    print *,'the value of zeta(s) at the last iteration of s:',Riemann_zeta(s)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    print *, 'Test 12 : computing the non-trivial zero of zeta(s) near 0.5+i*401.8 with Newtons method'
    z=0.5_qp+401.83922860053321653991130437828723072_qp*ii
    s=0.5_qp+401.8_qp*ii  
    do k=1,6
        h=Riemann_zeta_prime(s)
        s=s-h(1)/h(2)
        print *,'next iteration of Newtons method:', s
    end do
    print *,'error:', abs(s-z)
    print *,'the value of zeta(s) at the last iteration of s:',Riemann_zeta(s)      
!#######################################################################
    end program test
