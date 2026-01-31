program bayesian_iteration
        implicit none
        integer, parameter :: n = 256, m = 961
        integer :: i, j, unit, ios, itr
        real(8) :: mean, rms, thresh, sumdph
        real(8) :: dph(n), mdph(n), shadows(m,n), image(m), noise(n), ratio(n), correction(m), shafac(m), D(m), maxD
        character(len=200) :: head, line
        real(8), parameter :: e = 0.001d0, tfac = 3.0d0

        !Opening and reading the detector plane histogram file.
        !----------------------------------------------------------
        open(unit=11, file="dph1.txt", action="read", status="old", iostat=ios)
        if (ios /= 0) then
                print*, "Error in opening the dph text file", ios
                stop
        end if

        read(11, *, iostat=ios) dph
        if (ios /= 0) then
                print*, "Error in reading the dph text file", ios
                stop
        end if
        close(11)
        !---------------------------------------------------------------

        !Opening and reading the model image (intensity) file.
        !----------------------------------------------------
        open(unit=13, file="image_file_1.txt", action="read", status="old", iostat=ios) !the text file is an array of 961 elements; all 1s representing initial image
        if (ios /= 0) then
                print*, "Error in opening the image text file", ios
                stop
        end if

        read(13, *, iostat=ios) image
        if (ios /= 0) then
                print*, "Error in reading the image text file", ios
                stop
        end if
        close(13)
        image = image/sum(image)
        !---------------------------------------------------------------

        !Adding noise and then normalizing the dph.
        !----------------------------------------------------------
        call random_number(noise) 
        noise = 0.05d0*noise
        dph = dph + noise
        sumdph = sum(dph)
        dph = dph/sum(dph)
        !----------------------------------------------------------

        !Reading the shadow arrays from its text file, normalizing and storing them in shadows(961x256) matrix.
        !------------------------------------------------------------------------------------------------------- 
        open(unit=14, file='shadow_arrays.txt', status='old', action='read') 
        open(unit=15, file='shadows_matrix.txt', status='replace', action='write')
        do i = 1, 3
                read(14, *) head
        end do
        do i = 1, 961
                read(14, *) line
                read(14, *) shadows(i,:)
                write(15, '(256(F2.0,1X))') shadows(i, :)
                shadows(i,:) = shadows(i,:)/sum(shadows(i,:))
                shafac(i) = 1.0d0/sum(shadows(i,:))
        end do
        close(14)
        close(15)
        !-----------------------------------------------------------------------------------------
        
        !Bayesian Iteration.
        !-----------------------------------------------------------------------------------------
        maxD = 5
        itr = 0
        do while (maxD > e .and. itr<1000)
                itr = itr + 1
                mdph = 0.0d0
                do i = 1, 961
                        mdph = mdph + image(i)*shadows(i,:)
                end do
                call random_number(noise)
                mdph = mdph + (0.01d0/256.d0)*noise
                mdph = mdph/sum(mdph)
                ratio = dph/mdph
                do i = 1, 961
                        correction(i) = sum(ratio*shadows(i,:))
                        image(i) = correction(i)*image(i)
                end do
                image = image/sum(image)

                mean = sum(image)/m
                rms = sqrt(sum((image-mean)**2)/m)
                thresh = tfac*rms + mean
                do i = 1, 961
                        if (image(i) > thresh) then
                                D(i) = abs(1 - correction(i))
                        else 
                                D(i) = 0.0d0
                        endif
                end do
        print*, "iteration no.:",itr,"max|1-Ci|:",maxval(D),"at",maxloc(D)
        maxD = maxval(D)
        end do
        !--------------------------------------------------------------------------

        !write(*, '(961(F9.6,1X))') correction
        print*, "Max. value of |1-Ci| is:", maxval(D)
        print*, "It occurs at:", maxloc(D)
        print*, "The total number of iterations is",itr
        open(unit=16, file='sky_profile.txt', status='replace', action='write')
        do i = 1, 961
                write(16,*) i, image(i)
        end do
        close(16)
        open(unit=17, file="final_image.txt", status="replace", action='write')
        image = (image*sumdph)*shafac
        do j = 1, 961
                write(17,'(I3,F11.7,1X)') j, image(j)
        end do
        close(17) 

end program bayesian_iteration

