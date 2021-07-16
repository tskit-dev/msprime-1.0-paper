all: ms discoal argon seq-gen simbac

ms:
	make -C msdir
	cp msdir/ms ./

seq-gen:
	rm -fR seq-gen Seq-Gen* 1.3.4.tar.gz
	wget https://github.com/rambaut/Seq-Gen/archive/1.3.4.tar.gz
	tar -zxf 1.3.4.tar.gz
	make -C Seq-Gen-1.3.4/source
	cp Seq-Gen-1.3.4/source/seq-gen ./

discoal:
	wget https://github.com/kern-lab/discoal/archive/v0.1.4.tar.gz
	tar -xf v0.1.4.tar.gz
	rm v0.1.4.tar.gz
	cd discoal-0.1.4/ && make
	cp discoal-0.1.4/discoal ./


simbac:
	git clone https://github.com/tbrown91/SimBac.git
	cd SimBac && g++ *.cpp -lgsl -lgslcblas -lm -O2 -o SimBac
	cp SimBac/SimBac ./simbac

argon:
	git clone https://github.com/pierpal/ARGON.git
	cp ARGON/ARGON.0.1.jar ./ARGON.jar

clean:
	rm -fR ms discoal* ARGON* msms* Seq-Gen* seq-gen
