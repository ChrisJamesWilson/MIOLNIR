c------------------------------------------------------------------------c
      program isolf_split
c------------------------------------------------------------------------c
c     Aaron Dotter 2006 - send comments to aaron.dotter@gmail.com        c
c *** Revised 2009 to handle much wider files (more columns)             c
c *** Revised 2008 to correct problem with filenames on some systems     c
c     If filenames need to be rounded upward, e.g. a10099 -> 10100       c
c     set fudge >0 in the parameter statement below.    (Thanks TAS)     c
c *** Revised 2007 to allow for younger ages, output file names differ   c
c     from the original version: age string is now 5 characters long     c
c --- This program takes an isochrone or luminosity function file        c
c     containing multiple ages and divides it up into multiple files     c
c     each containing one age.                                           c
c     To compile:    F77 -o isolf_split isolf_split.f                    c
c     It deals with almost everything as strings, so no numbers will     c
c     change in the process. The names of the new files will begin with  c
c     "a" for age, followed by five integers, followed by the root       c
c     name which is the name of the original file.                       c
c     examples:         15  Gyr from "test.cmd" = "a15000test.cmd"       c
c                       6.5 Gyr from "test.cmd" = "a06500test.cmd"       c
c                       250 Myr from "test.cmd" = "a00250test.cmd"       c
c     On systems that support the commands IArgC and getarg (g77),       c
c     the program may be run entirely from the command line as           c
c            ./isolf_split <filename>                                    c
c     by uncommenting the code under "run in quiet mode".                c
c------------------------------------------------------------------------c
      implicit none
      integer i,j,nages,max,in,out,iage
      parameter(max=100,in=1,out=2)
      integer eeps(max)
      character agestr*5,eepstr*6,nagestr*19,hdr1*53,hdr2*53,hdr3*60
      character dataline*700,hdr4*700,sep*60,filename*25,newfile*31
      character cage*5,imf*15
      double precision age(max), fudge
      parameter(fudge=0.0d0)
      logical lquiet
      data lquiet/.false./

c run in quiet mode
      if(IArgC().gt.0) then
         lquiet=.true.
         call getarg(1,filename)
      endif
     
c run interactively
      if(.not.lquiet)then
         print *, '-------------------------------------------'
         print *, '     Isochrone Age Separation Program'
         print *, '-------------------------------------------'
         print *, "  Enter filename:"
         read *, filename
         print *, '-------------------------------------------'
      endif

c output large iso or LF file containing >1 ages
      open(in,file=filename,status='old')
      read(in,'(a16,i2)') nagestr,nages
      read(in,'(a)') sep
      read(in,'(a)') hdr1
      read(in,'(a)') hdr2
      read(in,*)
      read(in,'(a)') hdr3
      read(in,*)
c loop through individual ages
      do i=1,nages
         read(in,'(a5,f6.3,a6,i3,a15)')agestr,age(i),eepstr,eeps(i),imf
         read(in,'(a)') hdr4
c create new filename by converting integer age to a 5-character string
         iage=int(age(i)*1.0d3+fudge)
         newfile(1:1)='a'
         if(iage.lt.10000.and.iage.ge.1000)then
            write(cage,'(i1,i4)') 0,iage
         else if(iage.lt.1000)then
            write(cage,'(i1,i1,i3)') 0,0,iage
         else
            write(cage,'(i5)') iage
         endif
         newfile(2:6)=cage
         newfile(7:)=filename
c write new isochrone to file
         open(out,file=newfile,status='unknown')
         write(out,'(a16,i2)') nagestr,1
         write(out,'(a)') sep
         write(out,'(a)') hdr1
         write(out,'(a)') hdr2
         write(out,'(a)') sep
         write(out,'(a)') hdr3
         write(out,'(a)') sep
         write(out,'(a5,f6.3,a6,i3)') agestr,age(i),eepstr,eeps(i)
         write(out,'(a)') hdr4(1:len_trim(hdr4))
         do j=1,eeps(i)
            read(in,'(a)') dataline
            write(out,'(a)') dataline(1:len_trim(dataline))
         enddo
         if(i.lt.nages) read(in,*)
         if(i.lt.nages) read(in,*)
         close(out)
      enddo
      close(in)

c all done - print result
      if(.not.lquiet)then
         print *, '-------------------------------------------'
         print *, '  Results output to a***',filename
         print *, '-------------------------------------------'
      endif

      end
c------------------------------------------------------------------------c
c------------------------------------------------------------------------c
