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
	wget https://github.com/kr-colab/discoal/archive/v0.1.5.tar.gz
	tar -xf v0.1.5.tar.gz
	rm v0.1.5.tar.gz
	sed -e 's/makeGametesMS/\/\/makeGametesMS/' discoal-0.1.5/discoal_multipop.c > tmp
	mv tmp discoal-0.1.5/discoal_multipop.c
	cd discoal-0.1.5/ &&  make 
	cp discoal-0.1.5/discoal ./


simbac:
	git clone https://github.com/tbrown91/SimBac.git
	cd SimBac && g++ *.cpp -lgsl -lgslcblas -lm -O2 -o SimBac
	cp SimBac/SimBac ./simbac

fastSimBac:
	git clone https://bitbucket.org/nicofmay/fastsimbac.git
	cd fastsimbac/fastSimBac_linux && make
	cp fastsimbac/fastSimBac_linux/fastSimBac ./

argon:
	git clone https://github.com/pierpal/ARGON.git
	cp ARGON/ARGON.0.1.jar ./ARGON.jar

clean:
	rm -fR ms discoal* ARGON* msms* Seq-Gen* seq-gen
