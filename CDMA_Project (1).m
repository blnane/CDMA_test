clear;
close all;

    %###Message1-BPSK Modulation###
%data = input('Enter Binary Data In Form of 0 and 1 in [ ] : ');
data1 = [1 0 0 1 1 0 1 1];
msg_signal = rectpulse(data1, 100);

data1(data1(:)==0) = -1; %NRZ form
msg_signal_NRZ1 = rectpulse(data1, 100);

N = length(data1);       %Message length
fc = 10;                %Carrier frequency
ac = 1;                 %Carrier amplitude
t1 = 0.01:0.01:N;       %Time period

bpsk_mod1 = sqrt(2/ac)*cos(2*pi*fc*t1);
bpsk_signal1 = msg_signal_NRZ1 .* bpsk_mod1;

    %###Message1 PN Generation, CDMA Modulation###
sr1 = [1 -1 1 -1];
pn1 = [];
for i = 1:N
    for j=1:10
        pn1 = [pn1 sr1(4)];

        if sr1(4) == sr1(3)
            temp = -1;
        else
            temp = 1;
        end

        sr1(4) = sr1(3);
        sr1(3) = sr1(2);
        sr1(2) = sr1(1);
        sr1(1) = temp;
    end
end

pn_upsampled1 = [];
for i = 1:length(pn1)
    for j = 1:10
        pn_upsampled1 = [pn_upsampled1 pn1(i)];
    end
end

sigtx1 = bpsk_signal1 .* pn_upsampled1;

    %###Message2-BPSK Modulation###
%data = input('Enter Binary Data In Form of 0 and 1 in [ ] : ');
data2 = [1 1 0 0 1 0 1 0];
msg_signal = rectpulse(data2, 100);

data2(data2(:)==0) = -1; %NRZ form
msg_signal_NRZ2 = rectpulse(data2, 100);

N = length(data2);      %Message length
fc = 10;                %Carrier frequency
ac = 1;                 %Carrier amplitude
t1 = 0.01:0.01:N;       %Time period

bpsk_mod2 = sqrt(2/ac)*cos(2*pi*fc*t1);
bpsk_signal2 = msg_signal_NRZ2 .* bpsk_mod2;

    %###Message2 PN Generation, CDMA Modulation###
sr2 = [-1 1 -1 1];
pn2 = [];
for i = 1:N
    for j=1:10
        pn2 = [pn2 sr2(4)];

        if sr2(4) == sr2(3)
            temp = -1;
        else
            temp = 1;
        end

        sr2(4) = sr2(3);
        sr2(3) = sr2(2);
        sr2(2) = sr2(1);
        sr2(1) = temp;
    end
end

pn_upsampled2 = [];
for i = 1:length(pn2)
    for j = 1:10
        pn_upsampled2 = [pn_upsampled2 pn2(i)];
    end
end

sigtx2 = bpsk_signal2 .* pn_upsampled2;

    %###CDMA Message Combine###
sigtx = sigtx1 + sigtx2;

    %###Make noise###
sigtonoise = 10;
composite_signal_10 = awgn(sigtx, sigtonoise);
sigtonoise = 5;
composite_signal_5 = awgn(sigtx, sigtonoise);
sigtonoise = 0;
composite_signal_0 = awgn(sigtx, sigtonoise);

    %###Message 1 CDMA, BPSK Demodulation###
sigtonoise = 20;
composite_signal = awgn(sigtx, sigtonoise);

rs1 = composite_signal .* pn_upsampled1; %CDMA Demodulation

rm1 = []; %BPSK Demodulation
bpsk_demod = rs1 .* bpsk_mod1;
for i = 1:100:size(bpsk_mod1, 2)
    rm1 = [rm1 trapz(t1(i:i+99), bpsk_demod(i:i+99))];
end
rm1(rm1(:)<=0) = -1;
rm1(rm1(:)>=0) = 1;

rm_signal1 = rectpulse(rm1, 100);

    %###Message 2 CDMA, BPSK Demodulation###
sigtonoise = 20;
composite_signal = awgn(sigtx, sigtonoise);

rs2 = composite_signal .* pn_upsampled2; %CDMA Demodulation

rm2 = []; %BPSK Demodulation
bpsk_demod = rs2 .* bpsk_mod2;
for i = 1:100:size(bpsk_mod2, 2)
    rm2 = [rm2 trapz(t1(i:i+99), bpsk_demod(i:i+99))];
end
rm2(rm2(:)<=0) = -1;
rm2(rm2(:)>=0) = 1;

rm_signal2 = rectpulse(rm2, 100);

%=========OUTPUT==========
    %###Message-BPSK Modulation output###
figure('Name', 'Message BPSK Modulation', 'NumberTitle', 'off');

subplot(2,2,1); plot(msg_signal_NRZ1);
axis([0 length(msg_signal_NRZ1) -1.2 1.2]);
title('Message 1 Signal in NRZ form');
xlabel('n');
ylabel('x(n)');
grid on;

subplot(2,2,2); plot(msg_signal_NRZ2);
axis([0 length(msg_signal_NRZ2) -1.2 1.2]);
title('Message 2 Signal in NRZ form');
xlabel('n');
ylabel('x(n)');
grid on;

subplot(2,2,3); plot(bpsk_signal1);
axis([0 100*N -2 2]);
title('Message 1 BPSK signal');
xlabel('Time Period(t)');
ylabel('x(t)');
grid on;

subplot(2,2,4); plot(bpsk_signal2);
axis([0 100*N -2 2]);
title('Message 2 BPSK signal');
xlabel('Time Period(t)');
ylabel('x(t)');
grid on;

%###PN Generation, CDMA Modulation output###
figure('Name', 'PN Generation and CDMA', 'NumberTitle', 'off');

subplot(3,2,1); stem(pn_upsampled1);
axis([0 length(pn_upsampled1) -1.2 1.2]);
title(sprintf('Message 1\nPN sequence for data Upsampled'));
xlabel('n');
ylabel('x(n)');
grid on;

subplot(3,2,2); stem(pn_upsampled2);
axis([0 length(pn_upsampled2) -1.2 1.2]);
title(sprintf('Message 2\nPN sequence for data Upsampled'));
xlabel('n');
ylabel('x(n)');
grid on;

subplot(3,2,3); plot(sigtx1);
title('Message 1 CDMA Signal');
xlabel('Time Period(t)');
ylabel('x(t)');
grid on;

subplot(3,2,4); plot(sigtx2);
title('Message 2 CDMA Signal');
xlabel('Time Period(t)');
ylabel('x(t)');
grid on;

subplot(3,2,[5,6]); plot(sigtx);
title('Combined CDMA Signal');
xlabel('Time Period(t)');
ylabel('x(t)');
grid on;

%###Make Noise output###
figure('Name', 'Different SNR', 'NumberTitle', 'off');

subplot(3,1,1); plot(composite_signal_10);
title('SNR = 10db');
xlabel('Time Period(t)');
ylabel('x(t)');
grid on;

subplot(3,1,2); plot(composite_signal_5);
title('SNR = 5db');
xlabel('Time Period(t)');
ylabel('x(t)');
grid on;

subplot(3,1,3); plot(composite_signal_0);
title('SNR = 0db');
xlabel('Time Period(t)');
ylabel('x(t)');
grid on;

%###Message 1 CDMA, BPSK Demodulation output###
figure('Name', 'Message 1 CDMA Reciever', 'NumberTitle', 'off');

subplot(2,2,1); plot(composite_signal);
title(sprintf('Tx signal + noise\n SNR=20db'));
xlabel('Time Period(t)');
ylabel('x(t)');
grid on;

subplot(3,2,2); stem(pn_upsampled1);
axis([0 length(pn_upsampled1) -1.2 1.2]);
title(sprintf('Message 1\nPN sequence for data Upsampled'));
xlabel('n');
ylabel('x(n)');
grid on;

subplot(2,2,3); plot(rs1);
title('CDMA Demodulated signal');
xlabel('Time Period(t)');
ylabel('x(t)');
grid on;

subplot(2,2,4); plot(rm_signal1);
axis([0 length(rm_signal1) -1.2 1.2]);
title('Recieved Message signal in NRZ');
xlabel('n');
ylabel('x(t)');
grid on;

%###Message 2 CDMA, BPSK Demodulation output###
figure('Name', 'Message 2 CDMA Reciever', 'NumberTitle', 'off');

subplot(2,2,1); plot(composite_signal);
title(sprintf('Tx signal + noise\n SNR=20db'));
xlabel('Time Period(t)');
ylabel('x(t)');
grid on;

subplot(3,2,2); stem(pn_upsampled2);
axis([0 length(pn_upsampled2) -1.2 1.2]);
title(sprintf('Message 2\nPN sequence for data Upsampled'));
xlabel('n');
ylabel('x(n)');
grid on;

subplot(2,2,3); plot(rs2);
title('CDMA Demodulated signal');
xlabel('Time Period(t)');
ylabel('x(t)');
grid on;

subplot(2,2,4); plot(rm_signal2);
axis([0 length(rm_signal2) -1.2 1.2]);
title('Recieved Message signal in NRZ');
xlabel('n');
ylabel('x(t)');
grid on;