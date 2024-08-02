%Javier Arribas 2015
% Flexiband I 1 binary read file
clear all;
close all;
flag_initial_packet=1;
fileID = fopen('C:\signals\test_ifen_I_1b\cap3_I_1b.bin');

sample_index=1;

%Flexiband I 1b
%             <FrontendMode name="I-1b (L1_E1abc_B1)">
%             <Variant major="1" minor="1" micro="1"/>
%             <bitfile name="flexiband_I-1b.bit"/>
%             <usbAltIntfc>1</usbAltIntfc>
%             <framing preamble="0x55,0xAA" counter="4" crc="4"/>
%             <vgas>
%                 <vga2>L1_E1abc</vga2>
%             </vgas>
%             <FrontendBand name="L1_E1abc" bits="4" samples="2" complex="1" usb="1" vga="2" filterIdx="0" sampleRate="40"/>
%             </FrontendMode>

header_offset_bytes=6;
crc_bytes=4;
%while (~feof(fileID))
for (n=1:1:100)
    disp(num2str(n));
    [packet,count] = fread(fileID,1024,'uint8');
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
            for (n_sample=1:1:(1024-header_offset_bytes-crc_bytes))
                sample_byte=packet(header_offset_bytes+n_sample);
                Q=int8(bitand(sample_byte,15,'uint8'));
                I=int8(bitshift(sample_byte,-4,'uint8'));
                if I>=8
                    I=I-16;
                end
                if Q>=8
                    Q=Q-16;
                end
                I=2*I+1;
                Q=2*Q+1;
                s(sample_index)=double(I)+1i*double(Q);
                sample_index=sample_index+1;
            end
%             crc_hex_str=['0x' dec2hex(uint8(packet(n_sample+header_offset_bytes+1))) ...
%                 ' 0x' dec2hex(uint8(packet(n_sample+header_offset_bytes+2))) ...
%                 ' 0x' dec2hex(uint8(packet(n_sample+header_offset_bytes+3))) ...
%                 ' 0x' dec2hex(uint8(packet(n_sample+header_offset_bytes+4)))]
        else
             disp(['Header incorrect!, received this header ' header_hex_str]);
        end
    end
end
disp('End of file reached');
fclose(fileID);

%analysis

%s=exp(i*2*pi*(1/Fs)*205e3*(1:1:1000));
Fs=40e6;

pwelch(s,hamming(1024),[],[],Fs,'centered');

  figure
  y=read_complex_char_binary('C:\signals\test_ifen_I_1b\cap3_I_1b_L1_E1abc.bin',10000);
  plot(cumsum(abs(y-s(1:10000).')))
% [a,f]=pwelch(s,[],[],[],Fs,'twosided');
% f_twosides=[(-Fs+f(1+length(f)/2:end))' f(1:length(f)/2)'];
% plot(f_twosides,[a(1+length(f)/2:end)' a(1:length(f)/2)']);
% 
% figure
% 
% semilogy(f_twosides,[a(1+length(f)/2:end)' a(1:length(f)/2)']);
% 
% 
% s_cut=s(1:1014*3*8*40);
% figure;
% plot(real(s_cut));
% hold on;
% plot(imag(s_cut),'r');

%figure;

%plot(atan2(real(s_cut),imag(s_cut)));

%figure;
%plot(atan2(imag(s_cut),real(s_cut)));

% figure;
% plot(cumsum(atan2(imag(s_cut),real(s_cut))));

