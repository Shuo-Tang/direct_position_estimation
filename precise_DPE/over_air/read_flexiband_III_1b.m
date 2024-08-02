%Javier Arribas 2015
% Flexiband I 1 binary read file
clear all;
close all;
flag_initial_packet=1;
fileID = fopen('C:\signals\test_ifen_III_1b\cap3_III_1b.bin');

sample_index1=1;
sample_index2=1;
sample_index3=1;

%     <FrontendMode name="III-1b (L1/L2/L5 extended)">
%         <Variant major="3" minor="1" micro="1"/>
%         <bitfile name="flexiband_III-1b.bit"/>
%         <usbAltIntfc>1</usbAltIntfc>
%         <framing preamble="0x55,0xAA" counter="4" crc="6"/>
%         <vgas>
%             <vga1>L2</vga1>
%             <vga2>L1</vga2>
%             <vga3>L5</vga3>
%         </vgas>
%         <FrontendBand name="L2" bits="4" samples="2" complex="1" usb="1" vga="1" filterIdx="0" sampleRate="20"/>
%         <FrontendBand name="L1" bits="4" samples="2" complex="1" usb="1" vga="2" filterIdx="1" sampleRate="20"/>
%         <FrontendBand name="L5" bits="4" samples="4" complex="1" usb="1" vga="3" filterIdx="2" sampleRate="40"/>
%     </FrontendMode>

header_offset_bytes=6;
crc_bytes=6;
n_samples_per_frame=252; %(1024-6-6)/4
sample_mux_bytes=4;
enable_data_decoding=1;
%while (~feof(fileID))
for n=1:1:100
    [packet,count] = fread(fileID,1024,'uint8');
    disp(num2str(n));
    if count>0 
        header_hex_str=['0x' dec2hex(uint8(packet(1))) ' 0x' dec2hex(uint8(packet(2)))];
        if (packet(1)==85 && packet(2)==170)
           counter=packet(3)*2^24+packet(4)*2^16+packet(5)*2^8+packet(6);
           if (flag_initial_packet==1)
               flag_initial_packet=0;
           else
               if (last_counter~=counter-1)
                   disp(['Missing packet: last_counter=' num2str(last_counter) ' and counter=' num2str(counter)]);
               end
           end
           last_counter=counter;
           %decode packet
           if (enable_data_decoding==1)
                for n_sample=0:1:(n_samples_per_frame)

                    %L1
                    sample_byte=packet(header_offset_bytes+n_sample*sample_mux_bytes+1);
                    Q2=int8(bitand(sample_byte,15,'uint8'));
                    I2=int8(bitshift(sample_byte,-4,'uint8'));
                    if I2>=8
                        I2=I2-16;
                    end
                    if Q2>=8
                        Q2=Q2-16;
                    end
                    Q2=2*Q2+1;
                    I2=2*I2+1;
                    s2(sample_index2)=double(I2)+1i*double(Q2);
                    sample_index2=sample_index2+1;

                    %L2
                    sample_byte=packet(header_offset_bytes+n_sample*sample_mux_bytes+2);
                    Q1=int8(bitand(sample_byte,15,'uint8'));
                    I1=int8(bitshift(sample_byte,-4,'uint8'));
                    if I1>=8
                        I1=I1-16;
                    end
                    if Q1>=8
                        Q1=Q1-16;
                    end
                    I1=2*I1+1;
                    Q1=2*Q1+1;
                    s1(sample_index1)=double(I1)+1i*double(Q1);
                    sample_index1=sample_index1+1;

                    %L5
                    sample_byte=packet(header_offset_bytes+n_sample*sample_mux_bytes+3);
                    Q3=int8(bitand(sample_byte,15,'uint8'));
                    I3=int8(bitshift(sample_byte,-4,'uint8'));
                    if I3>=8
                        I3=I3-16;
                    end
                    if Q3>=8
                        Q3=Q3-16;
                    end
                    Q3=2*Q3+1;
                    I3=2*I3+1;
                    s3(sample_index3)=double(I3)+1i*double(Q3);
                    sample_index3=sample_index3+1;      

                    sample_byte=packet(header_offset_bytes+n_sample*sample_mux_bytes+4);
                    Q3=int8(bitand(sample_byte,15,'uint8'));
                    I3=int8(bitshift(sample_byte,-4,'uint8'));
                    if I3>=8
                        I3=I3-16;
                    end
                    if Q3>=8
                        Q3=Q3-16;
                    end
                    I3=2*I3+1;
                    Q3=2*Q3+1;
                    s3(sample_index3)=double(I3)+1i*double(Q3);
                    sample_index3=sample_index3+1;  
                end
%                  crc_hex_str=['0x' dec2hex(uint8(packet(1024-crc_bytes+1))) ...
%                       ' 0x' dec2hex(uint8(packet(1024-crc_bytes+2))) ...
%                       ' 0x' dec2hex(uint8(packet(1024-crc_bytes+3))) ...
%                       ' 0x' dec2hex(uint8(packet(1024-crc_bytes+4))) ...
%                       ' 0x' dec2hex(uint8(packet(1024-crc_bytes+5))) ...
%                       ' 0x' dec2hex(uint8(packet(1024-crc_bytes+6)))]
               end
        else
             disp(['Header incorrect!, received this header ' header_hex_str]);
        end
    end
end
disp('End of file reached');
fclose(fileID);



%write_complex_binary(s1,'c:\signals\cap_L1a.dat');
%write_complex_binary(s2,'c:\signals\cap_L2a.dat');
%write_complex_binary(s3,'c:\signals\cap_L5a.dat');

% %analysis
% 

y1=read_complex_char_binary('C:\signals\test_ifen_III_1b\cap3_III_1b__L1.bin',1000);
plot(cumsum(abs(y1-s1(1:1000).')))
figure
  Fs1=20e6;
  pwelch(s1,hamming(1024),[],[],Fs1,'centered');
%  [a1,f1]=pwelch(s1,[],[],[],Fs1,'twosided');
%  f_twosides1=[(-Fs1+f1(1+length(f1)/2:end))' f1(1:length(f1)/2)'];
%  plot(f_twosides1,[a1(1+length(f1)/2:end)' a1(1:length(f1)/2)']);
%  
  figure
%  

y2=read_complex_char_binary('C:\signals\test_ifen_III_1b\cap3_III_1b__L2.bin',1000);
plot(cumsum(abs(y2-s2(1:1000).')))
 Fs2=20e6;
 figure
  pwelch(s2,hamming(1024),[],[],Fs2,'centered');
%  f_twosides2=[(-Fs2+f2(1+length(f2)/2:end))' f2(1:length(f2)/2)'];
%  plot(f_twosides2,[a2(1+length(f2)/2:end)' a2(1:length(f2)/2)']);
%  
  figure
%  
   Fs3=40e6;
  y3=read_complex_char_binary('C:\signals\test_ifen_III_1b\cap3_III_1b__L5.bin',1000);
  pwelch(s3,hamming(1024),[],[],Fs3,'centered');
  figure
  plot(cumsum(abs(y3-s3(1:1000).')))
%  f_twosides3=[(-Fs3+f3(1+length(f3)/2:end))' f3(1:length(f3)/2)'];
%  plot(f_twosides3,[a3(1+length(f3)/2:end)' a3(1:length(f3)/2)']);
%  
%  figure
% % 
%  semilogy(f_twosides1,[a1(1+length(f1)/2:end)' a1(1:length(f1)/2)']);
%  
%   figure
% % 
%  semilogy(f_twosides2,[a2(1+length(f2)/2:end)' a2(1:length(f2)/2)']);
%  
%   figure
% % 
%  semilogy(f_twosides3,[a3(1+length(f3)/2:end)' a3(1:length(f3)/2)']);
% 
% 
%  s1_cut=s1(1:1014);
%  figure;
%  plot(real(s1_cut));
%  hold on;
%  plot(imag(s1_cut),'r');
% 
% %figure;
% 
% %plot(atan2(real(s_cut),imag(s_cut)));
% 
% %figure;
% %plot(atan2(imag(s_cut),real(s_cut)));
% 
% figure;
% plot(cumsum(atan2(imag(s_cut),real(s_cut))));

