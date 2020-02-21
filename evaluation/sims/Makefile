all: ms msms scrm discoal argon

ms:
	make -C msdir
	cp msdir/ms ./

msms:
	wget https://www.mabs.at/ewing/msms/msms3.2rc-b163.zip
	unzip msms3.2rc-b163.zip
	cp msms/lib/msms.jar ./

scrm:
	wget https://github.com/scrm/scrm/releases/download/v1.7.2/scrm-src.tar.gz
	tar -xf scrm-src.tar.gz
	cd scrm-1.7.2 && ./configure && make
	cp scrm-1.7.2/scrm ./

discoal:
	wget https://github.com/kern-lab/discoal/archive/v0.1.4.tar.gz
	tar -xf v0.1.4.tar.gz
	cd discoal-0.1.4/ && make 
	cp discoal-0.1.4/discoal ./

argon:
	git clone https://github.com/pierpal/ARGON.git
	cp ARGON/ARGON.0.1.jar ./ARGON.jar

clean:
	rm -fR ms scrm* discoal* v0* ARGON* msms*
