module collision_nuclear
    use calc_df 
    implicit none
    
    private

    ! Collision type modelled
    LOGICAL, PARAMETER   :: do_dd_fusion = .true.
    
    REAL(num), PARAMETER :: eps = EPSILON(1.0_num)
    REAL(num), PARAMETER :: one_m_2eps = 1.0_num - 2.0_num * eps
    REAL(num), PARAMETER :: one_p_2eps = 1.0_num + 2.0_num * eps
    
    REAL(num), PARAMETER :: e_rest = m0 * c**2
    REAL(num), PARAMETER :: e_rest_ev = e_rest / ev
    REAL(num), PARAMETER :: cc = c**2

    REAL(num), PARAMETER :: pi4_eps2_c4 = 4.0_num * pi * epsilon0**2 * c**4
    REAL(num), PARAMETER :: two_thirds = 2.0_num / 3.0_num
    REAL(num), PARAMETER :: pi_fac = &
         (4.0_num * pi / 3.0_num)**(1.0_num / 3.0_num)
    REAL(num), PARAMETER :: j_to_kev = 6.241509e+15 
    
    public binary_rkine
    type binary_rkine
       type(particle_list) :: p_list  
       real(num)           :: mass, charge, weight
       real(num)           :: user_factor
       real(num)           :: dens
       real(num)           :: log_lambda
       integer             :: n_species
       integer             :: he_species
     contains
       procedure         :: initialize
       procedure         :: do_collide
       procedure         :: model_dd_fusion
       procedure         :: create_fusion_products
       procedure         :: get_species_id       
       procedure, nopass :: generate_produced_particle
    END TYPE binary_rkine

  CONTAINS
    !> initialisation procedure
    SUBROUTINE initialize(this)
      CLASS(binary_rkine), INTENT(INOUT)    ::   this
      INTEGER :: ispecies

      DO ispecies = 1, n_species
         IF (species_list(ispecies)%name == TRIM('neutron')) THEN
                this%n_species = ispecies
             ELSE IF  (species_list(ispecies)%name == TRIM('helium')) THEN
                this%he_species = ispecies
          END IF  
       END DO
       !PRINT*,'nuclear_reaction species id, n: ', this%n_species, ' he: ', this%he_species 
      
     END SUBROUTINE initialize

     SUBROUTINE get_species_id(this, name, id)
       CLASS(binary_rkine), INTENT(INOUT)    :: this
       CHARACTER(len = * )                   :: name
       INTEGER, INTENT(OUT)                  :: id
       INTEGER                               :: ispecies
       
       DO ispecies = 1, n_species
          IF (species_list(ispecies)%name == TRIM(name)) THEN
             id = ispecies
          END IF
       END DO

     END SUBROUTINE get_species_id

    SUBROUTINE do_collide(this, n_list, he_list)
      CLASS(binary_rkine), INTENT(INOUT)    :: this
      TYPE(particle_list), INTENT(INOUT)    :: n_list, he_list      

      IF ( do_dd_fusion ) THEN
         PRINT *, "n_reactions -> do_dd_fusion: species mass: ", this%mass
         call this%model_dd_fusion(n_list, he_list)
      ELSE
         PRINT *, "n_reactions -> no reactions selected ... "
      ENDIF
    END SUBROUTINE do_collide
     

    SUBROUTINE model_dd_fusion(this, n_list, he_list)
      CLASS(binary_rkine), INTENT(INOUT)    :: this
      TYPE(particle),      POINTER          :: current, impact
      TYPE(particle_list), INTENT(INOUT)      :: n_list, he_list      
      TYPE(particle),      POINTER          :: n1_part, n2_part, he1_part, he2_part
      REAL(num) :: factor
      INTEGER(i8) :: icount, k, pcount, npairs, n2max, n_ratio
      REAL(num) :: ran1, ran2, s12, cosp, sinp, s_fac, v_rel
      REAL(num) :: sinp_cos, sinp_sin, s_prime, s_fac_prime
      REAL(num) :: a, a_inv, p_perp, p_tot, v_sq, gamma_rel_inv
      REAL(num) :: p_perp2, p_perp_inv, cell_fac
      REAL(num), DIMENSION(3) :: p1, p2, p3, p4, vc, v1, v2, p5, p6,v12
      REAL(num), DIMENSION(3) :: p1_norm, p2_norm
      REAL(num), DIMENSION(3,3) :: mat
      REAL(num) :: p_mag, p_mag2, fac, gc, vc_sq
      REAL(num) :: gm1, gm2, gm3, gm4, gm, gc_m1_vc
      REAL(num) :: m1, m2, q1, q2, m3, m4, dens_23
      REAL(num) :: minW, maxW, wp, cell_v, prob_fusion    
      REAL(num) :: v12_norm, ek_r , cs 
      LOGICAL   :: e_ok
      
      ! Case of pairing with one specie (so called intra collisions) 
      ! Number of m_particles in the list
      icount = this%p_list%count
      IF (icount <= 1) RETURN      
      ! Number of collisions
      pcount = icount / 2 + MOD(icount, 2_i8)
      !PRINT *, "model_dd_fusion -> : processing ncollisions: ", pcount, " non-even excess: ", MOD(icount, 2_i8)
      
      ! Use per-species particle properties 
      m1 = this%mass 
      m2 = this%mass 
      q1 = this%charge
      q2 = this%charge
      ! NRatio factor ( Higginson eq 17 for like-particles )
      ! NRatio_(like_particles) = N_pairs/Sampled_Pairs = (N*(N-1)/2)/(N/2) = N-1 
      n_ratio = icount - 1
      ! Cell volume (2D) units m^2
      cell_v = dx*dy
      
      !PRINT *, "model_dd_fusion ->  q1:", q1, " m1: ", m1
      !PRINT *, "model_dd_fusion -> NRatio: ", n_ratio, "cell_v:", cell_v

      ! Join tail to the head of the list to make it circular
      this%p_list%tail%next => this%p_list%head

      ! Initialise current and impact particle
      ! consecutively in the part. list
      current => this%p_list%head
      impact => current%next
      
      ! Loop over pairs   
      DO k = 1, pcount      
         ! Get Min and Max Weight
         minW= min(current%weight,impact%weight)
         maxW= max(current%weight,impact%weight)
         ! Momentums for paired particles
         p1 = current%part_p 
         p2 = impact%part_p 
         p1_norm = p1 / m0
         p2_norm = p2 / m0
         ! Remove non-moving particles 
         IF (DOT_PRODUCT(p1_norm, p1_norm) < eps &
              .AND. DOT_PRODUCT(p2_norm, p2_norm) < eps) CYCLE         
         ! Remove particles with the same momentum
         vc = (p1_norm - p2_norm)
         IF (DOT_PRODUCT(vc, vc) < eps) CYCLE         
         p1_norm = p1 / m1
         p2_norm = p2 / m2
         v12 = p1_norm - p2_norm
         v12_norm = SQRT(DOT_PRODUCT(v12, v12))  

         ! Ek_r = 1/2 * m12 * v12**2
         ek_r = 0.5_num * m1 * 0.5_num * (v12_norm**2)
         ek_r = ek_r * j_to_kev         
         cs = xsec( ek_r )

         ! Probability for fusion
         prob_fusion = (n_ratio * maxW * v12_norm * cs * dt) / cell_v
         prob_fusion = prob_fusion * this%user_factor 
         !PRINT*,"model_dd_fusion  prob_fusion: ", prob_fusion, "ek_r: ", ek_r,  "xsec: ",  cs, ' dt_step: ', dt 
         ! Adjusting probability in seldom cases where P_fusion > 1 
         DO WHILE (prob_fusion>.99_num)       
            this%user_factor = this%user_factor/10.0_num
            prob_fusion = (this%user_factor * (n_ratio * maxW * v12_norm * cs * dt)) / cell_v
            PRINT*," Adjusting prob_fusion: ", prob_fusion   
         END DO
           
         ! Check if fusion event occurs
         ran1 = random()
         IF ( ran1 < prob_fusion ) THEN
            !PRINT*,' Fusion event has occured ran1: ', ran1, ' Pfusion: ', prob_fusion
            ! Create fusion products
            CALL this%create_fusion_products(p1, m1, p2, m2, p3, m3, p4, m4, e_ok)
            IF (e_ok) THEN
               !PRINT*, ' created Fusion products m1 ', m1, ' pxyz:', p1(1), p1(2), p1(3)
               !PRINT*, ' created Fusion products m2 ', m2, ' pxyz:', p2(1), p2(2), p2(3)              
               !PRINT*, ' created Fusion products m3 ', m3, ' pxyz:', p3(1), p3(2), p3(3)
               !PRINT*, ' created Fusion products m4 ', m4, ' pxyz:', p4(1), p4(2), p4(3)
            END IF

          ! Calculate new weight for fusion products 
          wp = minW / this%user_factor 
          ! Reduce weight of both reactant   
          current%weight = current%weight - wp
          impact%weight  = impact%weight  - wp

          ! Creates product particles
          ! Here we take by defaut the position of the test particle (current)
          ! In principle it would be better to correct the weigth by a factor
          ! ( q1,2/q1+q2) and make the particle in both particle p1, p2 positions
          ! But for that we need to be able to assign the  total charge to the
          ! created charged products ( in this case He)
          CALL generate_produced_particle(current%part_pos, p3, wp, n1_part)
          !CALL generate_produced_particle(impact%part_pos, p3, 0.5_num*wp, n2_part)          
          CALL generate_produced_particle(current%part_pos, p4, wp, he1_part)
          !CALL generate_produced_particle(impact%part_pos, p4, 0.5_num*wp, he2_part)                              

          ! adding particles to corresponding lists 
          CALL add_particle_to_partlist(n_list, n1_part)
          !CALL add_particle_to_partlist(n_list, n2_part)
          CALL add_particle_to_partlist(he_list, he1_part)
          !CALL add_particle_to_partlist(he_list, he2_part)          
          
          NULLIFY(n1_part)
          !NULLIFY(n2_part)
          NULLIFY(he1_part)
          !NULLIFY(he2_part)

         END IF
         
         ! Get the next part. pair 
         current => impact%next
         impact => current%next
#ifdef PREFETCH
         CALL prefetch_particle(current)
         CALL prefetch_particle(impact)
#endif
      END DO ! do k, npairs  
      

      ! restore the tail of the list
      NULLIFY(this%p_list%tail%next)      
      
    END SUBROUTINE model_dd_fusion


    SUBROUTINE create_fusion_products(this, p1, m1, p2, m2, p3, mass3, p4, mass4, conserv)
      CLASS(binary_rkine), INTENT(INOUT)    :: this
      REAL(num), DIMENSION(3), INTENT(IN)   :: p1, p2
      REAL(num)              , INTENT(IN)   :: m1, m2
      REAL(num), DIMENSION(3), INTENT(OUT)  :: p3, p4
      REAL(num),               INTENT(OUT)  :: mass3, mass4 
      REAL(num), DIMENSION(3)               :: va, vb, u, V, Vcm, Vprime, a_vec
      REAL(num), DIMENSION(3)               :: p1_norm, p2_norm, p3_norm, p4_norm 
      REAL(num), DIMENSION(3)               :: V3prime, V4prime, vv3, vv4, u3, u4, v3, v4      
      REAL(num), DIMENSION(3,3)             :: mat, t_mat
      REAL(NUM)                             :: cosphi, sinphi
      REAL(NUM)                             :: costheta, sintheta
      REAL(NUM)                             :: mr, V3mag, vab
      REAL(NUM)                             :: cosA, cosP, sinA, sinP, phi      
      REAL(NUM)                             :: gm1, gm2, gm3, gm4
      LOGICAL  ,                INTENT(OUT) :: conserv            
      REAL(NUM), PARAMETER                  :: Q  = 3.269e3 !(keV)
      REAL(NUM), PARAMETER                  :: m3 = 1838.7 * m0
      REAL(NUM), PARAMETER                  :: m4 = 5497.9 * m0


      ! Fusion Reaction Kinematics
      ! Following algorithm described in 
      ! D. P. Higginson, A. Link, A. Schmidt, J. Comput. Phys. 388, 439 (2019)

      
      ! Pre-collision velocities
      va = p1 / m1
      vb = p2 / m2

      ! Step A1
      u = va - vb
      vab = SQRT(DOT_PRODUCT(u,u))
      
      ! Step A2
      IF (SQRT(u(1)**2 + u(2)**2 ) > 0) THEN
         cosphi = u(1) / SQRT(u(1)**2 + u(2)**2 ) 
         sinphi = u(2) / SQRT(u(1)**2 + u(2)**2 ) 
      ELSE
         cosphi = 0.0_num
         sinphi = 1.0_num
      ENDIF
      
      costheta = u(3) / SQRT(DOT_PRODUCT(u,u))
      sintheta = SQRT((u(1)**2 + u(2)**2) / DOT_PRODUCT(u,u))
      
      mat(1,1) =  cosphi*costheta 
      mat(1,2) =  sinphi*costheta 
      mat(1,3) = -1.*sintheta  
      mat(2,1) = -1.*sinphi    
      mat(2,2) =  cosphi    
      mat(2,3) =  0.0_num 
      mat(3,1) =  cosphi*sintheta 
      mat(3,2) =  sinphi*sintheta 
      mat(3,3) =  costheta

      V = MATMUL(mat, u)
      
      ! Step A3
      Vcm = V * (m1/(m1+m2))
      Vprime = V - Vcm

      ! Step B
      mr = (m1*m2)/(m1+m2)      
      ! DB Check me
      V3mag = SQRT( ( 2._num/(m3+m4) ) * (m4/m3) * ( 0.5_num*mr*((vab)**2) + (Q/j_to_kev)) )
      !PRINT*, 'V3mag= ', V3mag 
      
      cosA = 1._num - 2._num  * random() 
      sinA = SQRT( 1._num - cosA**2) 
      phi = 2._num * pi * random()
      cosP = COS(phi)
      sinP = SIN(phi)
      
      a_vec(1) = sinA*cosP
      a_vec(2) = sinA*sinP
      a_vec(3) = cosA
      
      V3prime = V3mag * a_vec 
      V4prime = -1.* SQRT(m3/m4) * V3prime
      
      ! Step C1
      vv3 = V3prime + Vcm
      vv4 = V4prime + Vcm
      
      ! Step C2
      t_mat = TRANSPOSE(mat)
      u3 = MATMUL(t_mat, vv3)
      u4 = MATMUL(t_mat, vv4)
      
      ! Step C3
      v3 = u3 + vb
      v4 = u4 + vb
      
      mass3 = m3
      mass4 = m4
      p3 = mass3 * v3 ! SI units (kg. m. s-1)
      p4 = mass4 * v4 

      ! Check conservation of total Energy
      ! using  gamma = sqrt (1 + ( p_vec/m/c)^2)
      p1_norm = p1 / m1 / c
      p2_norm = p2 / m2 / c
      p3_norm = p3 / m3 / c
      p4_norm = p4 / m4 / c
      
      gm1 = SQRT(DOT_PRODUCT(p1_norm, p1_norm) + 1.0_num) * m1
      gm2 = SQRT(DOT_PRODUCT(p2_norm, p2_norm) + 1.0_num) * m2
      gm3 = SQRT(DOT_PRODUCT(p3_norm, p3_norm) + 1.0_num) * m3
      gm4 = SQRT(DOT_PRODUCT(p4_norm, p4_norm) + 1.0_num) * m4
      
      !PRINT*, ' Total energy conservation: gm1+gm2: ' , (gm1+gm2)*c**2 , ' gm3+gm4: ', (gm3+gm4)*c**2 
      IF ( ABS( (gm1+gm2)*c**2 - (gm3+gm4)*c**2 ) < 1.e-8 ) THEN
       conserv = .true.
      END IF  
       
      
    END SUBROUTINE create_fusion_products

    FUNCTION xsec(e)
      ! Bosh-Hale fusion cross section
      IMPLICIT NONE
      REAL(num) :: xsec
      REAL(num), INTENT(IN) :: e ! keV
      REAL(num), PARAMETER :: BG = 31.3970 ! sqrt(keV)
      REAL(num), PARAMETER :: A1 =  5.3701e4
      REAL(num), PARAMETER :: A2 =  3.3027e2
      REAL(num), PARAMETER :: A3 = -1.2706e-1
      REAL(num), PARAMETER :: A4 =  2.9327e-5
      REAL(num), PARAMETER :: A5 = -2.5151e-9
      REAL(num), PARAMETER :: B1 =  0.0
      REAL(num), PARAMETER :: B2 =  0.0
      REAL(num), PARAMETER :: B3 =  0.0
      REAL(num), PARAMETER :: B4 =  0.0
      REAL(num)            :: se, se_top, se_down

      ! Bosh-Hale parametrisation
      !         sigma(E) = S(E)/(E*exp(BG/sqrt(E)))

      IF ( e > 0.01) THEN ! threshold dE > 10 eV
         se_top = A1 + e*(A2 + e*(A3 + e*(A4 + e*A5))) 
         se_down = 1. + e*(B1 + e*(B2 + e*(B3 + e*B4)))            
         se = se_top / se_down
         xsec = (se / (e * exp( BG / SQRT(e) ) ) )*1e-31 ! m^2
      ELSE
         xsec = 0.0_num
      END IF

    END FUNCTION xsec

    SUBROUTINE generate_produced_particle(pos, p, weight, produced_part)
      REAL(num), DIMENSION(3), INTENT(IN)   :: p
      REAL(num), DIMENSION(c_ndims), INTENT(IN) :: pos
      REAL(num), INTENT(IN) :: weight
      TYPE(particle), POINTER, INTENT(OUT) :: produced_part
      
      CALL create_particle(produced_part)
      produced_part%part_pos = pos
      produced_part%part_p = p
      produced_part%weight = weight
      
    END SUBROUTINE generate_produced_particle
    
    SUBROUTINE model_dd_fusion2(p_list, mass, charge, weight, &
         dens, log_lambda, user_factor)
      ! Perform collisions between particles of the same species.
      TYPE(particle_list), INTENT(IN) :: p_list
      REAL(num), INTENT(IN) :: mass, charge, weight
      REAL(num), INTENT(IN) :: user_factor
      REAL(num), INTENT(IN) :: dens, log_lambda
      TYPE(particle), POINTER :: current, impact
      REAL(num) :: factor
      INTEGER(i8) :: icount, k, pcount
      REAL(num) :: ran1, ran2, s12, cosp, sinp, s_fac, v_rel
      REAL(num) :: sinp_cos, sinp_sin, s_prime, s_fac_prime
      REAL(num) :: a, a_inv, p_perp, p_tot, v_sq, gamma_rel_inv
      REAL(num) :: p_perp2, p_perp_inv, cell_fac
      REAL(num), DIMENSION(3) :: p1, p2, p3, p4, vc, v1, v2, p5, p6
      REAL(num), DIMENSION(3) :: p1_norm, p2_norm
      REAL(num), DIMENSION(3,3) :: mat
      REAL(num) :: p_mag, p_mag2, fac, gc, vc_sq
      REAL(num) :: gm1, gm2, gm3, gm4, gm, gc_m1_vc
      REAL(num) :: m1, m2, q1, q2, dens_23
      REAL(num), PARAMETER :: pi4_eps2_c4 = 4.0_num * pi * epsilon0**2 * c**4
      REAL(num), PARAMETER :: two_thirds = 2.0_num / 3.0_num
      REAL(num), PARAMETER :: pi_fac = &
           (4.0_num * pi / 3.0_num)**(1.0_num / 3.0_num)
      factor = 0.0_num
      
      ! Intra-species collisions
      icount = p_list%count
      
      ! If there aren't enough particles to collide, then don't bother
      IF (icount <= 1) RETURN
      
      ! Number of collisions
      pcount = icount / 2 + MOD(icount, 2_i8)
      
#ifdef PER_SPECIES_WEIGHT
      ! Factor of 2 due to intra species collisions
      ! See Section 4.1 of Nanbu
      factor = user_factor / (pcount * weight * 2.0_num)
#else
      ! temporarily join tail to the head of the list to make it circular
      p_list%tail%next => p_list%head
      
      current => p_list%head
      impact => current%next
      DO k = 1, pcount
         factor = factor + MIN(current%weight, impact%weight)
         current => impact%next
         impact => current%next
#ifdef PREFETCH
         CALL prefetch_particle(current)
         CALL prefetch_particle(impact)
#endif
      END DO
      factor = user_factor / factor / 2.0_num
#endif
      ! If possible, use per-species properties
      m1 = mass
      m2 = mass
      q1 = charge
      q2 = charge
      
      current => p_list%head
      impact => current%next
      
      ! Per-cell constant factors
      cell_fac = dens**2 * dt * n_coll_steps * factor * dx * dy 
      s_fac = cell_fac * log_lambda / pi4_eps2_c4
      dens_23 = dens**two_thirds
      s_fac_prime = cell_fac * pi_fac / dens_23
      
      ! Set some constants before looping on pairs
      !nuclear_reaction%na_23 = dens_23 
      !nuclear_reaction%nb_23 = dens_23
      
    
    DO k = 1, pcount
     ! Store kinetics of particle pair  
      !nuclear_reaction%pa  = current
      !nuclear_reaction%pb  = impact
      
#ifdef PER_PARTICLE_CHARGE_MASS
      m1 = current%mass
      m2 = impact%mass
      q1 = current%charge
      q2 = impact%charge
      !nuclear_reaction%ma  = current%mass
      !nuclear_reaction%mb  = impact%mass      
      !nuclear_reaction%mab = current%mass / impact%mass                  
#endif
      !nuclear_reaction%w_min = min(current%weight,impact%weight)
      !nuclear_reaction%w_max = max(current%weight,impact%weight)      
      ! DB: checkme : one need to correct time step due to diff. in weight.
      !nuclear_reaction%dt_corr = dt  
            
      p1 = current%part_p / c
      p2 = impact%part_p / c

      p1_norm = p1 / m0
      p2_norm = p2 / m0

      ! Two stationary particles can't collide, so don't try
      IF (DOT_PRODUCT(p1_norm, p1_norm) < eps &
          .AND. DOT_PRODUCT(p2_norm, p2_norm) < eps) CYCLE

      ! Ditto for two particles with the same momentum
      vc = (p1_norm - p2_norm)
      IF (DOT_PRODUCT(vc, vc) < eps) CYCLE

      p1_norm = p1 / m1
      gm1 = SQRT(DOT_PRODUCT(p1_norm, p1_norm) + 1.0_num) * m1

      p2_norm = p2 / m2
      gm2 = SQRT(DOT_PRODUCT(p2_norm, p2_norm) + 1.0_num) * m2

      gm = gm1 + gm2

      ! Pre-collision velocities
      v1 = p1 / gm1
      v2 = p2 / gm2

      ! Velocity of centre-of-mass (CM) expressed in the reference frame (Lab)
      vc = (p1 + p2) / gm
      vc_sq = DOT_PRODUCT(vc, vc)

      gamma_rel_inv = SQRT(1.0_num - vc_sq)
      gc = 1.0_num / gamma_rel_inv

      ! Store cm_v expressed in the Lab Frane
      !nuclear_reaction%cm_v = vc

      if (vc_sq < 1.e-6 ) then
         !nuclear_reaction%cm_g = 1.0_num + 0.5_num * vc_sq
         !nuclear_reaction%par1 = 0.5
      else
         !nuclear_reaction%cm_g = gc
         !nuclear_reaction%par1= (nuclear_reaction%cm_g - 1.0_num)/ vc_sq
      end if

      ! Lorentz transformed momentum of test particle
      gc_m1_vc = (gc - 1.0_num) / vc_sq
      p3 = p1 + (gc_m1_vc * DOT_PRODUCT(vc, v1) - gc) * gm1 * vc

      ! Lorentz tranformed gamma_1* gamma_2*
      v_sq = DOT_PRODUCT(vc, v1)
      gm3 = (1.0_num - v_sq) * gc * gm1
      v_sq = DOT_PRODUCT(vc, v2)
      gm4 = (1.0_num - v_sq) * gc * gm2

      ! CM Momentum norm  p* = p1* = p2*
      p_mag2 = DOT_PRODUCT(p3, p3)
      p_mag = SQRT(p_mag2)
      
      fac = (q1 * q2)**2 * s_fac / (gm1 * gm2)
      s12 = fac * gc * p_mag * c / gm * (gm3 * gm4 / p_mag2 + 1.0_num)**2

      ! Cold plasma upper limit for s12
      v_rel = gm * p_mag * c / (gm3 * gm4 * gc)
      s_prime = s_fac_prime * (m1 + m2) * v_rel / MAX(m1, m2)

      s12 = MIN(s12, s_prime)

      ran1 = random()
      ran2 = random() * 2.0_num * pi

      ! Inversion from Perez et al. PHYSICS OF PLASMAS 19, 083104 (2012)
      IF (s12 < 0.1_num) THEN
        cosp = 1.0_num + s12 * LOG(MAX(ran1, 5e-9_num))
      ELSE IF (s12 >= 0.1_num .AND. s12 < 3.0_num) THEN
        a_inv = 0.0056958_num + (0.9560202_num + (-0.508139_num &
             + (0.47913906_num + (-0.12788975_num + 0.02389567_num &
             * s12) * s12) * s12) * s12) * s12
        a = 1.0_num / a_inv
        cosp = a_inv * LOG(EXP(-a) + 2.0_num * ran1 * SINH(a))
      ELSE IF (s12 >= 3.0_num .AND. s12 < 6.0_num) THEN
        a = 3.0_num * EXP(-s12)
        cosp = LOG(EXP(-a) + 2.0_num * ran1 * SINH(a)) / a
      ELSE
        cosp = 2.0_num * ran1 - 1.0_num
      END IF

      ! Branches 2 and 3 can result in rounding errors
      cosp = MAX(MIN(cosp, 1.0_num), -1.0_num)

      sinp = SIN(ACOS(cosp))

      ! Calculate new momenta according to rotation by angle p
      p_perp2 = p3(1)**2 + p3(2)**2
      p_perp = SQRT(p_perp2)
      p_tot = SQRT(p_perp2 + p3(3)**2)
      p_perp_inv = 1.0_num / (p_perp + c_tiny)

      mat(1,1) =  p3(1) * p3(3) * p_perp_inv
      mat(1,2) = -p3(2) * p_tot * p_perp_inv
      mat(1,3) =  p3(1)
      mat(2,1) =  p3(2) * p3(3) * p_perp_inv
      mat(2,2) =  p3(1) * p_tot * p_perp_inv
      mat(2,3) =  p3(2)
      mat(3,1) = -p_perp
      mat(3,2) =  0.0_num
      mat(3,3) =  p3(3)

      sinp_cos = sinp * COS(ran2)
      sinp_sin = sinp * SIN(ran2)

      p3(1) = mat(1,1) * sinp_cos + mat(1,2) * sinp_sin + mat(1,3) * cosp
      p3(2) = mat(2,1) * sinp_cos + mat(2,2) * sinp_sin + mat(2,3) * cosp
      p3(3) = mat(3,1) * sinp_cos + mat(3,2) * sinp_sin + mat(3,3) * cosp

      p4 = -p3

      p5 = (p3 + (gc_m1_vc * DOT_PRODUCT(vc, p3) + gm3 * gc) * vc) * c
      p6 = (p4 + (gc_m1_vc * DOT_PRODUCT(vc, p4) + gm4 * gc) * vc) * c

      
      ! Update particle properties
      current%part_p = p5
      impact%part_p = p6

      current => impact%next
      impact => current%next
#ifdef PREFETCH
      CALL prefetch_particle(current)
      CALL prefetch_particle(impact)
#endif
   END DO ! do k, pcount
   
    ! restore the tail of the list
    NULLIFY(p_list%tail%next)

  END SUBROUTINE model_dd_fusion2

    
END MODULE collision_nuclear
