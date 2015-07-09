#### Test result
Successfully built /apps/sources/EasyBuild/custom/ngs-utils-15.07.1.eb

#### Overview of tested easyconfigs (in order)
 * **SUCCESS** _ngs-utils-15.07.1.eb_ 

#### Time info
 * start: Wed, 08 Jul 2015 07:55:48 +0000 (UTC)
 * end: Wed, 08 Jul 2015 07:56:16 +0000 (UTC)

#### EasyBuild info
 * easybuild-framework version: 2.1.1
 * easybuild-easyblocks version: 2.1.1
 * command line:
```
eb -f ngs-utils-15.07.1.eb
```
 * full configuration (includes defaults):
```
--buildpath=/apps//.tmp/easybuild/builds/
--cleanup-builddir
--disable-allow-modules-tool-mismatch
--disable-avail-cfgfile-constants
--disable-avail-easyconfig-constants
--disable-avail-easyconfig-licenses
--disable-avail-easyconfig-templates
--disable-avail-module-naming-schemes
--disable-avail-modules-tools
--disable-avail-repositories
--disable-debug
--disable-dry-run
--disable-dry-run-short
--disable-experimental
--disable-hidden
--disable-ignore-osdeps
--disable-info
--disable-job
--disable-list-toolchains
--disable-logtostdout
--disable-module-only
--disable-pretend
--disable-quiet
--disable-recursive-module-unload
--disable-regtest
--disable-sequential
--disable-set-gid-bit
--disable-show-default-configfiles
--disable-show-default-moduleclasses
--disable-skip
--disable-skip-test-cases
--disable-sticky-bit
--disable-update-modules-tool-cache
--disable-upload-test-report
--external-modules-metadata=""
--force
--ignore-dirs=.git,.svn
--installpath=/apps/
--logfile-format=easybuild,easybuild-%(name)s-%(version)s-%(date)s.%(time)s.log
--module-naming-scheme=EasyBuildMNS
--module-syntax=Tcl
--moduleclasses="['base', 'bio', 'cae', 'chem', 'compiler', 'data', 'debugger', 'devel', 'geo', 'ide', 'lang', 'lib', 'math', 'mpi', 'numlib', 'perf', 'phys', 'system', 'toolchain', 'tools', 'vis']"
--modules-tool=Lmod
--repository=FileRepository
--repositorypath=/home/umcg-rkanninga/.local/easybuild/ebfiles_repo
--robot-paths=""
--sourcepath=/apps//sources/
--strict=warn
--subdir-modules=modules
--subdir-software=software
--suffix-modules-path=all
````

#### System info
 * _core count:_ 48
 * _cpu model:_ Intel(R) Xeon(R) CPU E5-2697 v2 @ 2.70GHz
 * _cpu speed:_ 2700.075
 * _cpu vendor:_ Intel
 * _gcc version:_ Using built-in specs.; Target: x86_64-redhat-linux; Configured with: ../configure --prefix=/usr --mandir=/usr/share/man --infodir=/usr/share/info --with-bugurl=http://bugzilla.redhat.com/bugzilla --enable-bootstrap --enable-shared --enable-threads=posix --enable-checking=release --with-system-zlib --enable-__cxa_atexit --disable-libunwind-exceptions --enable-gnu-unique-object --enable-languages=c,c++,objc,obj-c++,java,fortran,ada --enable-java-awt=gtk --disable-dssi --with-java-home=/usr/lib/jvm/java-1.5.0-gcj-1.5.0.0/jre --enable-libgcj-multifile --enable-java-maintainer-mode --with-ecj-jar=/usr/share/java/eclipse-ecj.jar --disable-libjava-multilib --with-ppl --with-cloog --with-tune=generic --with-arch_32=i686 --build=x86_64-redhat-linux; Thread model: posix; gcc version 4.4.7 20120313 (Red Hat 4.4.7-11) (GCC) ; 
 * _glibc version:_ 2.12
 * _hostname:_ calculon
 * _os name:_ SL
 * _os type:_ Linux
 * _os version:_ 6.6
 * _platform name:_ x86_64-unknown-linux
 * _python version:_ 2.6.6 (r266:84292, Jan 22 2014, 05:06:49) ; [GCC 4.4.7 20120313 (Red Hat 4.4.7-3)]
 * _system gcc path:_ /usr/bin/gcc
 * _system python path:_ /usr/bin/python

#### List of loaded modules
 * EasyBuild/2.1.1

#### Environment
```
BASH_ENV = /usr/share/lmod/lmod/init/bash
BASH_FUNC_module() = () {  eval $($LMOD_CMD bash "$@");
 [ $? = 0 ] && eval $(${LMOD_SETTARG_CMD:-:} -s sh)
}
CVS_RSH = ssh
EASYBUILD_BUILDPATH = /apps//.tmp/easybuild/builds/
EASYBUILD_INSTALLPATH = /apps/
EASYBUILD_MODULES_TOOL = Lmod
EASYBUILD_SOURCEPATH = /apps//sources/
EBDEVELEASYBUILD = /apps/software/EasyBuild/2.1.1/easybuild/EasyBuild-2.1.1-easybuild-devel
EBROOTEASYBUILD = /apps/software/EasyBuild/2.1.1
EBVERSIONEASYBUILD = 2.1.1
G_BROKEN_FILENAMES = 1
HISTCONTROL = ignoredups
HISTSIZE = 1000
HOME = /home/umcg-rkanninga
HOSTNAME = calculon
LANG = en_US.UTF-8
LD_LIBRARY_PATH = /apps/software/EasyBuild/2.1.1/lib
LESSOPEN = ||/usr/bin/lesspipe.sh %s
LIBRARY_PATH = /apps/software/EasyBuild/2.1.1/lib
LMOD_CMD = /usr/share/lmod/lmod/libexec/lmod
LMOD_COLORIZE = yes
LMOD_DEFAULT_MODULEPATH = /gcc/groups/gcc/home/rkanninga/modules:/apps/modules/bio:/apps/modules/all:/etc/modulefiles:/usr/share/modulefiles:/usr/share/Modules/modulefiles:/usr/share/modulefiles/Linux:/usr/share/modulefiles/Core:/usr/share/lmod/lmod/modulefiles/Core
LMOD_DIR = /usr/share/lmod/lmod/libexec
LMOD_FULL_SETTARG_SUPPORT = no
LMOD_PKG = /usr/share/lmod/lmod
LMOD_PREPEND_BLOCK = normal
LMOD_SETTARG_CMD = :
LMOD_arch = x86_64
LMOD_sys = Linux
LOADEDMODULES = EasyBuild/2.1.1
LOGNAME = umcg-rkanninga
LS_COLORS = rs=0:di=38;5;27:ln=38;5;51:mh=44;38;5;15:pi=40;38;5;11:so=38;5;13:do=38;5;5:bd=48;5;232;38;5;11:cd=48;5;232;38;5;3:or=48;5;232;38;5;9:mi=05;48;5;232;38;5;15:su=48;5;196;38;5;15:sg=48;5;11;38;5;16:ca=48;5;196;38;5;226:tw=48;5;10;38;5;16:ow=48;5;10;38;5;21:st=48;5;21;38;5;15:ex=38;5;34:*.tar=38;5;9:*.tgz=38;5;9:*.arj=38;5;9:*.taz=38;5;9:*.lzh=38;5;9:*.lzma=38;5;9:*.tlz=38;5;9:*.txz=38;5;9:*.zip=38;5;9:*.z=38;5;9:*.Z=38;5;9:*.dz=38;5;9:*.gz=38;5;9:*.lz=38;5;9:*.xz=38;5;9:*.bz2=38;5;9:*.tbz=38;5;9:*.tbz2=38;5;9:*.bz=38;5;9:*.tz=38;5;9:*.deb=38;5;9:*.rpm=38;5;9:*.jar=38;5;9:*.rar=38;5;9:*.ace=38;5;9:*.zoo=38;5;9:*.cpio=38;5;9:*.7z=38;5;9:*.rz=38;5;9:*.jpg=38;5;13:*.jpeg=38;5;13:*.gif=38;5;13:*.bmp=38;5;13:*.pbm=38;5;13:*.pgm=38;5;13:*.ppm=38;5;13:*.tga=38;5;13:*.xbm=38;5;13:*.xpm=38;5;13:*.tif=38;5;13:*.tiff=38;5;13:*.png=38;5;13:*.svg=38;5;13:*.svgz=38;5;13:*.mng=38;5;13:*.pcx=38;5;13:*.mov=38;5;13:*.mpg=38;5;13:*.mpeg=38;5;13:*.m2v=38;5;13:*.mkv=38;5;13:*.ogm=38;5;13:*.mp4=38;5;13:*.m4v=38;5;13:*.mp4v=38;5;13:*.vob=38;5;13:*.qt=38;5;13:*.nuv=38;5;13:*.wmv=38;5;13:*.asf=38;5;13:*.rm=38;5;13:*.rmvb=38;5;13:*.flc=38;5;13:*.avi=38;5;13:*.fli=38;5;13:*.flv=38;5;13:*.gl=38;5;13:*.dl=38;5;13:*.xcf=38;5;13:*.xwd=38;5;13:*.yuv=38;5;13:*.cgm=38;5;13:*.emf=38;5;13:*.axv=38;5;13:*.anx=38;5;13:*.ogv=38;5;13:*.ogx=38;5;13:*.aac=38;5;45:*.au=38;5;45:*.flac=38;5;45:*.mid=38;5;45:*.midi=38;5;45:*.mka=38;5;45:*.mp3=38;5;45:*.mpc=38;5;45:*.ogg=38;5;45:*.ra=38;5;45:*.wav=38;5;45:*.axa=38;5;45:*.oga=38;5;45:*.spx=38;5;45:*.xspf=38;5;45:
MAIL = /var/spool/mail/umcg-rkanninga
MANPATH = /usr/share/lmod/lmod/share/man
MODULEPATH = /gcc/groups/gcc/home/rkanninga/modules:/apps/modules/bio:/apps/modules/all:/etc/modulefiles:/usr/share/modulefiles:/usr/share/Modules/modulefiles:/usr/share/modulefiles/Linux:/usr/share/modulefiles/Core:/usr/share/lmod/lmod/modulefiles/Core
MODULEPATH_ROOT = /usr/share/modulefiles
MODULESHOME = /usr/share/lmod/lmod
PATH = /apps/modules:/apps/software/EasyBuild/2.1.1/bin:/usr/share/lmod/lmod/libexec:/usr/lib64/qt-3.3/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/home/umcg-rkanninga/bin
PREFIX_ROOT = /apps/
PWD = /apps/sources/EasyBuild/custom
PYTHONPATH = /apps/software/EasyBuild/2.1.1/lib/python2.6/site-packages
QTDIR = /usr/lib64/qt-3.3
QTINC = /usr/lib64/qt-3.3/include
QTLIB = /usr/lib64/qt-3.3/lib
ROAN_GCC_HOME = /gcc/groups/gcc/tmp01/rkanninga/
ROAN_HOME = /gcc/groups/gcc/home/rkanninga
SHELL = /bin/bash
SHLVL = 2
SSH_CLIENT = 195.169.22.212 33746 22
SSH_CONNECTION = 195.169.22.212 33746 195.169.22.247 22
SSH_TTY = /dev/pts/12
TERM = xterm-256color
TEST_EASYBUILD_MODULES_TOOL = Lmod
USER = umcg-rkanninga
_ = /usr/bin/python
_LMFILES_ = /apps/modules/all/EasyBuild/2.1.1
_ModuleTable001_ = X01vZHVsZVRhYmxlXz17WyJhY3RpdmVTaXplIl09MSxiYXNlTXBhdGhBPXsiL2djYy9ncm91cHMvZ2NjL2hvbWUvcmthbm5pbmdhL21vZHVsZXMiLCIvYXBwcy9tb2R1bGVzL2JpbyIsIi9hcHBzL21vZHVsZXMvYWxsIiwiL2V0Yy9tb2R1bGVmaWxlcyIsIi91c3Ivc2hhcmUvbW9kdWxlZmlsZXMiLCIvdXNyL3NoYXJlL01vZHVsZXMvbW9kdWxlZmlsZXMiLCIvdXNyL3NoYXJlL21vZHVsZWZpbGVzL0xpbnV4IiwiL3Vzci9zaGFyZS9tb2R1bGVmaWxlcy9Db3JlIiwiL3Vzci9zaGFyZS9sbW9kL2xtb2QvbW9kdWxlZmlsZXMvQ29yZSIsfSxbImNfcmVidWlsZFRpbWUiXT04NjQwMCxbImNfc2hvcnRUaW1lIl09ZmFsc2UsZmFtaWx5PXt9LGluYWN0aXZlPXt9LG1UPXtFYXN5QnVp
_ModuleTable002_ = bGQ9e1siRk4iXT0iL2FwcHMvbW9kdWxlcy9hbGwvRWFzeUJ1aWxkLzIuMS4xIixbImRlZmF1bHQiXT0wLFsiZnVsbE5hbWUiXT0iRWFzeUJ1aWxkLzIuMS4xIixbImxvYWRPcmRlciJdPTEscHJvcFQ9e30sWyJzaG9ydCJdPSJFYXN5QnVpbGQiLFsic3RhdHVzIl09ImFjdGl2ZSIsfSx9LG1wYXRoQT17Ii9nY2MvZ3JvdXBzL2djYy9ob21lL3JrYW5uaW5nYS9tb2R1bGVzIiwiL2FwcHMvbW9kdWxlcy9iaW8iLCIvYXBwcy9tb2R1bGVzL2FsbCIsIi9ldGMvbW9kdWxlZmlsZXMiLCIvdXNyL3NoYXJlL21vZHVsZWZpbGVzIiwiL3Vzci9zaGFyZS9Nb2R1bGVzL21vZHVsZWZpbGVzIiwiL3Vzci9zaGFyZS9tb2R1bGVmaWxlcy9MaW51eCIsIi91c3Ivc2hhcmUvbW9kdWxlZmlsZXMv
_ModuleTable003_ = Q29yZSIsIi91c3Ivc2hhcmUvbG1vZC9sbW9kL21vZHVsZWZpbGVzL0NvcmUiLH0sWyJzeXN0ZW1CYXNlTVBBVEgiXT0iL2V0Yy9tb2R1bGVmaWxlczovdXNyL3NoYXJlL21vZHVsZWZpbGVzOi91c3Ivc2hhcmUvTW9kdWxlcy9tb2R1bGVmaWxlczovdXNyL3NoYXJlL21vZHVsZWZpbGVzL0xpbnV4Oi91c3Ivc2hhcmUvbW9kdWxlZmlsZXMvQ29yZTovdXNyL3NoYXJlL2xtb2QvbG1vZC9tb2R1bGVmaWxlcy9Db3JlIixbInZlcnNpb24iXT0yLH0=
_ModuleTable_Sz_ = 3
```