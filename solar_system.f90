!!!**************************************************************
!!!*						                	                *
!!!*                     Solar System     	         	        *
!!!*							        	                    *
!!!*                Author: Dr. Andrea Boghi  		            *
!!!*							        	                    *
!!!*                       1: Sun                               *
!!!*                       2: Mercury                           *
!!!*                       3: Venus                             *
!!!*                       4: Earth                             *
!!!*                       5: Mars                              *
!!!*                       6: Jupyter                           *
!!!*                       7: Saturn                            *
!!!*                       8: Uranus                            *
!!!*                       9: Neptune                           *
!!!*                      10: Pluto                             *
!!!*                      11: Luna                              *
!!!*                      12: Ceres                             *
!!!*								                            *
!!!**************************************************************
    
Program solar_system
    implicit none
    integer n_b,nt
    logical go_forward
    parameter(n_b=12)
    real*8 M(n_b),x0(n_b),y0(n_b),z0(n_b),x(n_b),y(n_b),z(n_b)
    real*8 GM(n_b),u0(n_b),v0(n_b),w0(n_b),u(n_b),v(n_b),w(n_b)
    real*8 x_temp(n_b),y_temp(n_b),z_temp(n_b),u_temp(n_b),v_temp(n_b),w_temp(n_b)
    real*8 dxdt_1(n_b),dydt_1(n_b),dzdt_1(n_b),dudt_1(n_b),dvdt_1(n_b),dwdt_1(n_b)
    real*8 dxdt_2(n_b),dydt_2(n_b),dzdt_2(n_b),dudt_2(n_b),dvdt_2(n_b),dwdt_2(n_b)
    real*8 dxdt_3(n_b),dydt_3(n_b),dzdt_3(n_b),dudt_3(n_b),dvdt_3(n_b),dwdt_3(n_b)
    real*8 dxdt_4(n_b),dydt_4(n_b),dzdt_4(n_b),dudt_4(n_b),dvdt_4(n_b),dwdt_4(n_b)
    real*8 dt,GU,days_to_seconds, seconds_to_days, time, PI, rad_to_deg
    real*8 Omega(n_b),Omega0(n_b),Omega_temp(n_b), r(n_b), phi(n_b), theta(n_b)
    real*8 dOmegadt_1(n_b),dOmegadt_2(n_b),dOmegadt_3(n_b),dOmegadt_4(n_b),dOmegadt(n_b)
    real*8 drdt(n_b), dphidt(n_b), dthetadt(n_b)
    parameter(PI=4.d0*datan2(1.d0,1.d0))
    parameter(rad_to_deg=180.d0/PI)
    parameter(GU=6.672662873575023d-20) !!! [km**3] /( [kg] * [s**2] )
    parameter(days_to_seconds=60.d0 * 60.d0 * 24.d0) !!![s]
    parameter(seconds_to_days=1.d0/days_to_seconds) !!![days]
    CHARACTER(len=80)::f_name_0
    CHARACTER(len=70) :: cols_sun
    CHARACTER(len=47) :: cols_sc
    CHARACTER(len=67) :: cols_sv
    integer k,i,j

    open(1,file='imp_ss.dat')
    read(1,*) nt
    read(1,*) go_forward
    read(1,*) dt      !!! [days]
    read(1,*) f_name_0
    close(1)

    call define_mass(M,n_b)
    call define_initial_speed_and_position(x0,y0,z0,u0,v0,w0,f_name_0,n_b)

    dt = dt * days_to_seconds

    if(go_forward) then
        dt = dabs(dt)
    else
        dt = -dabs(dt)
    endif

    GM = GU * M

    !!! initialize cycle
    x = x0
    y = y0
    z = z0
    u = u0
    v = v0
    w = w0
    Omega = 0.d0

    cols_sun = '"t [days]","x [km]","y [km]","z [km]","u [km/s]","v [km/s]","w [km/s]"'
    cols_sc = ',"r [km]","lat [deg]","lon [deg]","Omega [rad]"'
    cols_sv = ',"drdt [km/s]","dlatdt [rad/s]","dlondt [rad/s]","dOmegadt [rad/s]"'

    !!! Time Loop
    open(11,file="Solar_orbit.csv")
    write(11,*) cols_sun

    open(12,file="Mercury_orbit.csv")
    write(12,*) cols_sun // cols_sc // cols_sv

    open(13,file="Venus_orbit.csv")
    write(13,*) cols_sun // cols_sc // cols_sv

    open(14,file="Earth_orbit.csv")
    write(14,*) cols_sun // cols_sc // cols_sv

    open(15,file="Mars_orbit.csv")
    write(15,*) cols_sun // cols_sc // cols_sv

    open(16,file="Jupiter_orbit.csv")
    write(16,*) cols_sun // cols_sc // cols_sv

    open(17,file="Saturn_orbit.csv")
    write(17,*) cols_sun // cols_sc // cols_sv

    open(18,file="Uranus_orbit.csv")
    write(18,*) cols_sun // cols_sc // cols_sv

    open(19,file="Neptune_orbit.csv")
    write(19,*) cols_sun // cols_sc // cols_sv

    open(20,file="Pluto_orbit.csv")
    write(20,*) cols_sun // cols_sc // cols_sv

    open(21,file="Earth_Moon_orbit.csv")
    write(21,*) cols_sun // cols_sc // cols_sv

    open(22,file="Ceres_orbit.csv")
    write(22,*) cols_sun // cols_sc // cols_sv

    do k=1,nt

        time = dt * (k -1)

        !!! Reset cycle
        x0 = x
        y0 = y
        z0 = z
        u0 = u
        v0 = v
        w0 = w
        Omega0 = Omega

        !!! First Step
        call DDT(n_b, GM, x0, y0, z0, u0, v0, w0, dxdt_1, dydt_1, dzdt_1, dudt_1, dvdt_1, dwdt_1, dOmegadt_1)

        u_temp = u0 + (dt/2.d0) * dudt_1
        v_temp = v0 + (dt/2.d0) * dvdt_1
        w_temp = w0 + (dt/2.d0) * dwdt_1
        x_temp = x0 + (dt/2.d0) * dxdt_1
        y_temp = y0 + (dt/2.d0) * dydt_1 
        z_temp = z0 + (dt/2.d0) * dzdt_1

        !!! Second Step
        call DDT(n_b, GM, x_temp, y_temp, z_temp, u_temp, v_temp, w_temp, dxdt_2, dydt_2, dzdt_2, dudt_2, dvdt_2, dwdt_2, &
                dOmegadt_2)

        u_temp = u0 + (dt/2.d0) * dudt_2
        v_temp = v0 + (dt/2.d0) * dvdt_2
        w_temp = w0 + (dt/2.d0) * dwdt_2
        x_temp = x0 + (dt/2.d0) * dxdt_2
        y_temp = y0 + (dt/2.d0) * dydt_2 
        z_temp = z0 + (dt/2.d0) * dzdt_2

        !!! Third Step
        call DDT(n_b, GM, x_temp, y_temp, z_temp, u_temp, v_temp, w_temp, dxdt_3, dydt_3, dzdt_3, dudt_3, dvdt_3, dwdt_3, &
                dOmegadt_3)

        u_temp = u0 + dt * dudt_3
        v_temp = v0 + dt * dvdt_3
        w_temp = w0 + dt * dwdt_3
        x_temp = x0 + dt * dxdt_3
        y_temp = y0 + dt * dydt_3 
        z_temp = z0 + dt * dzdt_3

        !!! Fourth Step
        call DDT(n_b, GM, x_temp, y_temp, z_temp, u_temp, v_temp, w_temp, dxdt_4, dydt_4, dzdt_4, dudt_4, dvdt_4, dwdt_4, &
                dOmegadt_4)

        u = u0 + ( dt / 6.d0 )*( dudt_1 + 2.d0 * dudt_2 + 2.d0 * dudt_3 + dudt_4 )
        v = v0 + ( dt / 6.d0 )*( dvdt_1 + 2.d0 * dvdt_2 + 2.d0 * dvdt_3 + dvdt_4 )
        w = w0 + ( dt / 6.d0 )*( dwdt_1 + 2.d0 * dwdt_2 + 2.d0 * dwdt_3 + dwdt_4 )
        x = x0 + ( dt / 6.d0 )*( dxdt_1 + 2.d0 * dxdt_2 + 2.d0 * dxdt_3 + dxdt_4 )
        y = y0 + ( dt / 6.d0 )*( dydt_1 + 2.d0 * dydt_2 + 2.d0 * dydt_3 + dydt_4 )
        z = z0 + ( dt / 6.d0 )*( dzdt_1 + 2.d0 * dzdt_2 + 2.d0 * dzdt_3 + dzdt_4 )
        Omega = Omega0 + ( dt / 6.d0 )*( dOmegadt_1 + 2.d0 * dOmegadt_2 + 2.d0 * dOmegadt_3 + dOmegadt_4 )

        call spherical_coordinates_and_velocities(n_b, x, y, z, u, v, w, r, theta, phi, drdt, dthetadt, dphidt, dOmegadt)

        write(11,101) time * seconds_to_days,',',x(1),',',y(1),',',z(1),',',u(1),',',v(1),',',w(1)
        
        write(12,102) time * seconds_to_days,',',x(2)-x(1),',',y(2)-y(1),',',z(2)-z(1),',',u(2)-u(1),',',v(2)-v(1),',', &
                        w(2)-w(1),',',r(2),',',rad_to_deg * theta(2),',',rad_to_deg * phi(2),',',Omega(2) &
                        ,',',drdt(2),',',dthetadt(2),',',dphidt(2),',',dOmegadt(2)
        
        write(13,102) time * seconds_to_days,',',x(3)-x(1),',',y(3)-y(1),',',z(3)-z(1),',',u(3)-u(1),',',v(3)-v(1),',', &
                        w(3)-w(1),',',r(3),',',rad_to_deg * theta(3),',',rad_to_deg * phi(3),',',Omega(3) &
                        ,',',drdt(3),',',dthetadt(3),',',dphidt(3),',',dOmegadt(3)
        
        write(14,102) time * seconds_to_days,',',x(4)-x(1),',',y(4)-y(1),',',z(4)-z(1),',',u(4)-u(1),',',v(4)-v(1),',', &
                        w(4)-w(1),',',r(4),',',rad_to_deg * theta(4),',',rad_to_deg * phi(4),',',Omega(4) &
                        ,',',drdt(4),',',dthetadt(4),',',dphidt(4),',',dOmegadt(4)
        
        write(15,102) time * seconds_to_days,',',x(5)-x(1),',',y(5)-y(1),',',z(5)-z(1),',',u(5)-u(1),',',v(5)-v(1),',', &
                        w(5)-w(1),',',r(5),',',rad_to_deg * theta(5),',',rad_to_deg * phi(5),',',Omega(5) &
                        ,',',drdt(5),',',dthetadt(5),',',dphidt(5),',',dOmegadt(5)
        
        write(16,102) time * seconds_to_days,',',x(6)-x(1),',',y(6)-y(1),',',z(6)-z(1),',',u(6)-u(1),',',v(6)-v(1),',', &
                        w(6)-w(1),',',r(6),',',rad_to_deg * theta(6),',',rad_to_deg * phi(6),',',Omega(6) &
                        ,',',drdt(6),',',dthetadt(6),',',dphidt(6),',',dOmegadt(6)

        write(17,102) time * seconds_to_days,',',x(7)-x(1),',',y(7)-y(1),',',z(7)-z(1),',',u(7)-u(1),',',v(7)-v(1),',', &
                        w(7)-w(1),',',r(7),',',rad_to_deg * theta(7),',',rad_to_deg * phi(7),',',Omega(7) &
                        ,',',drdt(7),',',dthetadt(7),',',dphidt(7),',',dOmegadt(7)

        write(18,102) time * seconds_to_days,',',x(8)-x(1),',',y(8)-y(1),',',z(8)-z(1),',',u(8)-u(1),',',v(8)-v(1),',', &
                        w(8)-w(1),',',r(8),',',rad_to_deg * theta(8),',',rad_to_deg * phi(8),',',Omega(8) &
                        ,',',drdt(8),',',dthetadt(8),',',dphidt(8),',',dOmegadt(8)

        write(19,102) time * seconds_to_days,',',x(9)-x(1),',',y(9)-y(1),',',z(9)-z(1),',',u(9)-u(1),',',v(9)-v(1),',', &
                        w(9)-w(1),',',r(9),',',rad_to_deg * theta(9),',',rad_to_deg * phi(9),',',Omega(9) &
                        ,',',drdt(9),',',dthetadt(9),',',dphidt(9),',',dOmegadt(9)

        write(20,102) time * seconds_to_days,',',x(10)-x(1),',',y(10)-y(1),',',z(10)-z(1),',',u(10)-u(1),',',v(10)-v(1),',', &
                        w(10)-w(1),',',r(10),',',rad_to_deg * theta(10),',',rad_to_deg * phi(10),',',Omega(10) &
                        ,',',drdt(10),',',dthetadt(10),',',dphidt(10),',',dOmegadt(10)

        write(21,102) time * seconds_to_days,',',x(11)-x(4),',',y(11)-y(4),',',z(11)-z(4),',',u(11)-u(4),',',v(11)-v(4),',', &
                        w(11)-w(4),',',r(11),',',rad_to_deg * theta(11),',',rad_to_deg * phi(11),',',Omega(11) &
                        ,',',drdt(11),',',dthetadt(11),',',dphidt(11),',',dOmegadt(11)

        write(22,102) time * seconds_to_days,',',x(12)-x(1),',',y(12)-y(1),',',z(12)-z(1),',',u(12)-u(1),',',v(12)-v(1),',', &
                        w(12)-w(1),',',r(12),',',rad_to_deg * theta(12),',',rad_to_deg * phi(12),',',Omega(12) &
                        ,',',drdt(12),',',dthetadt(12),',',dphidt(12),',',dOmegadt(12)

    enddo

    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    close(18)
    close(19)
    close(20)
    close(21)
    close(22)

    101 format(e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10)

    102 format(e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10,A,e18.10 &
                ,A,e18.10,A,e18.10,A,e18.10,A,e18.10)


end

subroutine spherical_coordinates_and_velocities(n_b, x, y, z, u, v, w, r, theta, phi, drdt, dthetadt, dphidt, dOmegadt)
    implicit none
    integer, INTENT(IN) :: n_b
    real*8, INTENT(IN) :: x(n_b), y(n_b), z(n_b), u(n_b), v(n_b), w(n_b)
    real*8, INTENT(OUT) :: r(n_b), theta(n_b), phi(n_b), drdt(n_b), dthetadt(n_b), dphidt(n_b), dOmegadt(n_b)
    integer i,j

    do i=1,n_b
        if(i==1) then
            drdt(i) = 0.d0
            dphidt(i) = 0.d0
            dthetadt(i) = 0.d0
            dOmegadt(i) = 0.d0
        else
            if(i==11) then
                j = 4
            else
                j = 1
            endif

            call spherical_coordinates(r(i), theta(i), phi(i), x(i), y(i), z(i), x(j), y(j), z(j))

            call spherical_velocities(r(i), phi(i), theta(i), u(i), v(i), w(i), u(j), v(j), w(j), &
                    drdt(i), dphidt(i), dthetadt(i), dOmegadt(i))

        endif

    enddo

    end

subroutine spherical_velocities(r, phi, theta, u, v, w, u0, v0, w0, drdt, dphidt, dthetadt, dOmegadt)
    implicit none
    real*8, INTENT(IN) :: r, phi, theta, u, v, w, u0, v0, w0
    real*8, INTENT(OUT) :: drdt, dphidt, dthetadt, dOmegadt
    real*8 e_r(3), e_phi(3), e_theta(3)

    e_r(1) = dcos(phi) * dcos(theta)
    e_r(2) = dsin(phi) * dcos(theta)
    e_r(3) = dsin(theta)
    e_phi(1) = -dsin(phi) 
    e_phi(2) = dcos(phi) 
    e_phi(3) = 0.d0
    e_theta(1) = -dcos(phi) * dsin(theta)
    e_theta(2) = -dsin(phi) * dsin(theta)
    e_theta(3) = dcos(theta)

    drdt = ( u - u0 ) * e_r(1) + ( v - v0 ) * e_r(2) + ( w - w0 ) * e_r(3)
    dphidt = ( ( u - u0 ) * e_phi(1) + ( v - v0 ) * e_phi(2) + ( w - w0 ) * e_phi(3) ) &
                        / ( r * dcos(theta) )
    dthetadt = ( ( u - u0 ) * e_theta(1) + ( v - v0 ) * e_theta(2) + ( w - w0 ) * &
                        e_theta(3) ) / r

    if(dphidt >= 0.d0) then
        dOmegadt = ( ( dthetadt )**2.d0 + ( dphidt * dcos( theta ) )**2.d0 )**0.5d0
    else
        dOmegadt = - ( ( dthetadt )**2.d0 + ( dphidt * dcos( theta ) )**2.d0 )**0.5d0
    endif

    end

subroutine spherical_coordinates(r, theta, phi, x, y, z, x0, y0, z0)
    implicit none
    real*8, intent(IN) :: x, x0, y, y0, z, z0
    real*8, intent(out) :: r, theta, phi
    real*8 r_xy

    call distance(r, x, x0, y, y0, z, z0)
    call distance(r_xy, x, x0, y, y0, 0.d0, 0.d0)

    phi = datan2( (y-y0) , (x-x0) )
    theta = datan2( (z-z0) , r_xy )

    end

subroutine distance(r, x, x0, y, y0, z, z0)
    implicit none
    real*8, intent(IN) :: x, x0, y, y0, z, z0
    real*8, intent(out) :: r

    r = ( (x - x0)**2.d0 + (y - y0)**2.d0 + (z - z0)**2.d0 )**0.5d0

    end

subroutine DDT(n_b, GM, x, y, z, u, v, w, dxdt, dydt, dzdt, dudt, dvdt, dwdt, dOmegadt)
    implicit none
    integer, INTENT(IN) :: n_b
    real*8, INTENT(IN) :: GM(n_b), x(n_b), y(n_b), z(n_b), u(n_b), v(n_b), w(n_b)
    real*8, INTENT(OUT) :: dxdt(n_b), dydt(n_b), dzdt(n_b), dudt(n_b), dvdt(n_b), dwdt(n_b), dOmegadt(n_b)
    integer i,j
    real*8 r(n_b), theta(n_b), phi(n_b), drdt(n_b), dphidt(n_b), dthetadt(n_b)
    real*8 a

    do i=1,n_b
        dudt(i) = 0.d0
        dvdt(i) = 0.d0
        dwdt(i) = 0.d0
        dxdt(i) = u(i)
        dydt(i) = v(i)
        dzdt(i) = w(i)
        do j=1,n_b

            if(i /= j) then

                call distance(r(i), x(i), x(j), y(i), y(j), z(i), z(j))

                a = - GM(j) / r(i)**2.d0

                dudt(i) = dudt(i) + a * ( x(i) - x(j) ) / r(i)
                dvdt(i) = dvdt(i) + a * ( y(i) - y(j) ) / r(i)
                dwdt(i) = dwdt(i) + a * ( z(i) - z(j) ) / r(i)

            endif

        enddo
    enddo

    call spherical_coordinates_and_velocities(n_b, x, y, z, u, v, w, r, theta, phi, drdt, dthetadt, dphidt, dOmegadt)

    end

subroutine define_mass(M,n_b)
    implicit none
    integer, INTENT(IN) :: n_b
    real*8, INTENT(OUT) :: M(n_b)

    !!! [kg]
    M(1) = 1.9885d+30
    M(2) = 3.302d+23
    M(3) = 4.8685d+24
    M(4) = 5.97219d+24
    M(5) = 6.4171d+23
    M(6) = 1.89818722d+27
    M(7) = 5.6834d+26
    M(8) = 8.6813d+25
    M(9) = 1.02409d+26
    M(10) = 1.307d+22
    M(11) = 7.349d+22
    M(12) = 9.1d+20

    end

subroutine define_initial_speed_and_position(x0,y0,z0,u0,v0,w0,f_name_0,n_b)
    implicit none
    integer, INTENT(IN) :: n_b
    real*8, INTENT(OUT) :: x0(n_b),y0(n_b),z0(n_b),u0(n_b),v0(n_b),w0(n_b)
    CHARACTER(len=80), intent(in) :: f_name_0

    open(1,file=f_name_0)
    read(1,*) x0(1)    !!! [km]
    read(1,*) y0(1)    !!! [km]
    read(1,*) z0(1)    !!! [km]
    read(1,*) u0(1)    !!! [km / s]
    read(1,*) v0(1)    !!! [km / s]
    read(1,*) w0(1)    !!! [km / s]
    read(1,*) x0(2)    !!! [km]
    read(1,*) y0(2)    !!! [km]
    read(1,*) z0(2)    !!! [km]
    read(1,*) u0(2)    !!! [km / s]
    read(1,*) v0(2)    !!! [km / s]
    read(1,*) w0(2)    !!! [km / s]
    read(1,*) x0(3)    !!! [km]
    read(1,*) y0(3)    !!! [km]
    read(1,*) z0(3)    !!! [km]
    read(1,*) u0(3)    !!! [km / s]
    read(1,*) v0(3)    !!! [km / s]
    read(1,*) w0(3)    !!! [km / s]
    read(1,*) x0(4)    !!! [km]
    read(1,*) y0(4)    !!! [km]
    read(1,*) z0(4)    !!! [km]
    read(1,*) u0(4)    !!! [km / s]
    read(1,*) v0(4)    !!! [km / s]
    read(1,*) w0(4)    !!! [km / s]
    read(1,*) x0(5)    !!! [km]
    read(1,*) y0(5)    !!! [km]
    read(1,*) z0(5)    !!! [km]
    read(1,*) u0(5)    !!! [km / s]
    read(1,*) v0(5)    !!! [km / s]
    read(1,*) w0(5)    !!! [km / s]
    read(1,*) x0(6)    !!! [km]
    read(1,*) y0(6)    !!! [km]
    read(1,*) z0(6)    !!! [km]
    read(1,*) u0(6)    !!! [km / s]
    read(1,*) v0(6)    !!! [km / s]
    read(1,*) w0(6)    !!! [km / s]
    read(1,*) x0(7)    !!! [km]
    read(1,*) y0(7)    !!! [km]
    read(1,*) z0(7)    !!! [km]
    read(1,*) u0(7)    !!! [km / s]
    read(1,*) v0(7)    !!! [km / s]
    read(1,*) w0(7)    !!! [km / s]
    read(1,*) x0(8)    !!! [km]
    read(1,*) y0(8)    !!! [km]
    read(1,*) z0(8)    !!! [km]
    read(1,*) u0(8)    !!! [km / s]
    read(1,*) v0(8)    !!! [km / s]
    read(1,*) w0(8)    !!! [km / s]
    read(1,*) x0(9)    !!! [km]
    read(1,*) y0(9)    !!! [km]
    read(1,*) z0(9)    !!! [km]
    read(1,*) u0(9)    !!! [km / s]
    read(1,*) v0(9)    !!! [km / s]
    read(1,*) w0(9)    !!! [km / s]
    read(1,*) x0(10)   !!! [km]
    read(1,*) y0(10)   !!! [km]
    read(1,*) z0(10)   !!! [km]
    read(1,*) u0(10)   !!! [km / s]
    read(1,*) v0(10)   !!! [km / s]
    read(1,*) w0(10)   !!! [km / s]
    read(1,*) x0(11)   !!! [km]
    read(1,*) y0(11)   !!! [km]
    read(1,*) z0(11)   !!! [km]
    read(1,*) u0(11)   !!! [km / s]
    read(1,*) v0(11)   !!! [km / s]
    read(1,*) w0(11)   !!! [km / s]
    read(1,*) x0(12)   !!! [km]
    read(1,*) y0(12)   !!! [km]
    read(1,*) z0(12)   !!! [km]
    read(1,*) u0(12)   !!! [km / s]
    read(1,*) v0(12)   !!! [km / s]
    read(1,*) w0(12)   !!! [km / s]
    close(1)

    !!! convert to [m]
    !x0 = x0 * 1.d+03
    !y0 = y0 * 1.d+03
    !z0 = z0 * 1.d+03

    !!! convert to [m / s]
    !u0 = u0 * 1.d+03
    !v0 = v0 * 1.d+03
    !w0 = w0 * 1.d+03

    end


