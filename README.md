# Signal-Generator

The signal generator is a mex-function for MATLAB that can be used to generate the response of a moving sound source and receiver in a reverberant environment. 
The user can specify the position of the source and the receiver at each discrete time instance. 
The output signal is computed by convolving the (anechoic) source signal with the time-varying room impulse response. 
Multiple receiver positions can be specified to generate multiple responses simultaneously. 
The room impulse responses are generated using the image method, proposed by Allen and Berkley in 1979 [1]. 
The user can control the reverberation time (or reflection coefficients), reflection order, room dimension and microphone directivity in a way similar to the RIR generator. 

This package includes a MATLAB example, the mex-function, and the source code of the mex-function.

More information can be found [here](https://www.audiolabs-erlangen.de/fau/professor/habets/software/signal-generator).

[1] J.B. Allen and D.A. Berkley, "Image method for efficiently simulating small-room acoustics," Journal Acoustic Society of America, 65(4), April 1979, p 943.
