      module integrate
        public :: initialize
        public :: calculateNext
        contains

        subroutine :: initialize
            !TODO
        end subroutine initialize

        subroutine calculateNext(x,y,z,vx,vy,vz,ax,ay,az,m,tm,n,i_energy,radius,cur_time)
            !cm - center mass, sq - squared, dst - distance
            integer :: n
            integer :: index_i          !first minimal distance body
            integer :: index_j          !second minimal distance body
            logical :: evaporated(n)    !flags that identify the evaporation status of a body
            real*8 :: radius            !cluster radius
            real*8 :: x(n)
            real*8 :: y(n)
            real*8 :: z(n)
            real*8 :: vx(n)
            real*8 :: vy(n)
            real*8 :: vz(n)
            real*8 :: ax(n)
            real*8 :: ay(n)
            real*8 :: az(n)
            real*8 :: m(n)
            real*8 :: evd_time(n)       !evaporated time
            real*8 :: tm                !total mass
            real*8 :: a_vector
            real*8 :: sq_min_dst
            real*8 :: sq_max_v
            real*8 :: dt
            real*8 :: i_energy          !initial energy
            real*8 :: k_energy          !kinetic energy
            real*8 :: p_energy          !potential energy
            real*8 :: m_energy          !merge energy
            real*8 :: e_coeff           !energy coefficient
            real*8 :: v_cm_x
            real*8 :: v_cm_y
            real*8 :: v_cm_z
            real*8 :: cm_x
            real*8 :: cm_y
            real*8 :: cm_z
            real*8 :: ctm               !center total mass
            real*8 :: cur_time          !current_time

            !finding maximum distance between two bodies
            sq_min_dst = 1d256

            do i = 1, n
                do j = i + 1, n
                    sq_min_dst = min(sq_min_dst,(x(i)-x(j))**2 +(y(i)-y(j))**2 +(z(i)-z(j))**2)
                    index_i = i
                    index_j = j
                end do
            end do

            !checking for merging

            if(sq_min_dst < 1d18) then

                !energy difference
                m_energy = -6.67408d-11*m(index_i)*m(index_j)/dsqrt(sq_min_dst)

                i_energy = i_energy - m_energy

                !new body velocity
                vx(index_i)=(m(index_i) * vx(index_i) + m(index_j) * vx(index_j))/&
                    (m(index_i) + m(index_j))
			    vy(index_i)=(m(index_i) * vy(index_i) + m(index_j) * vy(index_j))/&
                    (m(index_i) + m(index_j))
				vz(index_i)=(m(index_i) * vz(index_i) + m(index_j) * vz(index_j))/&
                    (m(index_i) + m(index_j))

                !new body position
                x(index_i) = (m(index_i)*x(index_i)+m(index_j)*x(index_j))/&
                    (m(index_i)+m(index_j))
				y(index_i) = (m(index_i)*y(index_i)+m(index_j)*y(index_j))/&
                    (m(index_i)+m(index_j))
				z(index_i) = (m(index_i)*z(index_i)+z(index_j)*z(index_j))/&
                    (m(index_i)+m(index_j))

                !new body mass
                m(index_i) = m(indexi) + m(index_j)

                do i = index_j, n-1
                    x(i)=x(i+1)
                    y(i)=y(i+1)
                    z(i)=z(i+1)
                    vx(i)=vx(i+1)
                    vy(i)=vy(i+1)
                    vz(i)=vz(i+1)
                end do

                n = n - 1

                return

            end if

            !finding the maximal velocity in the system
            sq_max_v = 0d0

            do i = 1, n
                sq_max_v = max(sq_max_v,vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
            end do

            !calculating time steps
            dt = dsqrt(sq_min_dst/sq_max_v)

            !calculating acceleration projections
            do i = 1, n
                ax(i) = 0d0
                ay(i) = 0d0
                az(i) = 0d0
                do j = 1, n
                    if(i.NE.j) then
                        a_vector = 6.67408d-11*m(j)/&
                            ((x(i)-x(j))**2 +(y(i)-y(j))**2 +(z(i)-z(j))**2)**(1.5d0)
                        ax(i) = ax(i) + a_vector*(x(j)-x(i))
                        ay(i) = ay(i) + a_vector*(y(j)-y(i))
                        az(i) = az(i) + a_vector*(z(j)-z(i))
                    end if
                end do
            end do

            !calculating new coordinates and velocity projections
            do i = 1, n
                x(i) = x(i)+vx(i)*dt+ax(i)*dt**2/2
                y(i) = y(i)+vy(i)*dt+ay(i)*dt**2/2
                z(i) = z(i)+vz(i)*dt+az(i)*dt**2/2
                vx(i) = vx(i)+ax(i)*dt
                vy(i) = vy(i)+ay(i)*dt
                vz(i) = vz(i)+az(i)*dt
            end do

            !updating current time

            cur_time = cur_time + dt

            !energy conservation

            k_energy = 0d0
            p_energy = 0d0

            !kinetic energy
            do i = 1, n
                k_energy = k_energy + 0.5d0*m(i)*(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
            end do

            !potential energy
            do i = 1, n
                    do j = i+1, n
                        p_energy = p_energy - &
                            6.67408d-11*m(i)*m(j)/&
                                dsqrt((x(i)-x(j))**2 +(y(i)-y(j))**2 +(z(i)-z(j))**2)
                    end do
            end do

            !energy coefficient
            e_coeff = dsqrt(1 + (i_energy-(k_energy+p_energy))/k_energy)

            do i = 1, n
                vx(i) = vx(i)*e_coeff
                vy(i) = vy(i)*e_coeff
                vz(i) = vz(i)*e_coeff
            end do

            !making the system static

            v_cm_x = 0d0
            v_cm_y = 0d0
            v_cm_z = 0d0

            !total impulse projections
            do i = 1, n
                v_cm_x = v_cm_x + m(i)*vx(i)
                v_cm_y = v_cm_y + m(i)*vy(i)
                v_cm_z = v_cm_z + m(i)*vz(i)
            end do

            !center mass velocity projections
            v_cm_x = v_cm_x/tm
            v_cm_y = v_cm_y/tm
            v_cm_z = v_cm_z/tm

            do i = 1, n
                vx(i) = vx(i) - v_cm_x
                vy(i) = vy(i) - v_cm_y
                vz(i) = vz(i) - v_cm_z
            end do

            !updating evaporation status

            cm_x = 0d0
            cm_y = 0d0
            cm_z = 0d0
            ctm = 0d0

            do i = 1, n
                if(.NOT.evaporated(i)) then
                    ctm = ctm + m(i)
                    cm_x = cm_x + m(i)*x(i)
                    cm_y = cm_y + m(i)*y(i)
                    cm_z = cm_z + m(i)*z(i)
                end if
            end do

            cm_x = cm_x/ctm
            cm_y = cm_y/ctm
            cm_z = cm_z/ctm

            do i = 1, n
                if(.NOT.evaporated(i)) then
                    if(dsqrt((x(i)-cm_x)**2 + (y(i)-cm_y)**2 + (z(i)-cm_z)**2).GT.3*radius) then
                        if(0.5d0*m(i)*(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))-6.67408d-11*m(i)&
                            *(ctm-m(i)/dsqrt((x(i)-cm_x)**2 + (y(i)-cm_y)**2 + (z(i)-cm_z)**2 )).GT.0) then
                                evaporated(i) = .TRUE.
                                evd_time(i) = cur_time
                        end if
                    end if
                end if
            end do


        end subroutine calculateNext
      end module integrate
