# Leaky-Reluplex
Leaky-Reluplex : The extended algorithm of Reluplex, which is used to verify the DNNs with Leaky ReLU activation functions.

#### The paper of Reluplex:
 https://arxiv.org/abs/1702.01135

#### The codes of Reluplex:
https://github.com/guykatzz/ReluplexCav2017

#### Compilation Instructions

The implementation was run and tested on Ubuntu 16.04.

 Compiling GLPK:

      cd glpk-4.60
      ./configure_glpk.sh
      make
      make install

 Compiling the Reluplex core:

      cd reluplex
      make

 Test the Leaky-Reluplex (keep in the folder /Leaky-Reluplex/reluplex ):

      ./test.sh

The test log can see the file ./test.txt

The test case in Leaky-Reluplex is the same as in Reluplex except that the k value of Leaky ReLU function is 0.2.