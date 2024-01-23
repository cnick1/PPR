function [A, B] = acrs_sod(C, M, N)
% Subroutine to perform SOD calculations
% This code is a MATLAB translation of the Fortran code ACRS from [1,2]. 
% References: [1] P. R. Graves-Morris and D. E. Roberts, "ACRSFORTRAN: a
%             subroutine and procedure for the rapid calculation of simple
%             off-diagonal rational approximants." Mendeley, Jan. 1975.
%             doi: 10.17632/TM9GZ2TGR6.1
%             [2] P. R. Graves-Morris and D. E. Roberts, "A subroutine and 
%             procedure for the rapid calculation of simple off-diagonal
%             rational approximants,” Computer Physics Communications, vol.
%             9, no. 1, pp. 46–50, Jan. 1975, 
%             doi: 10.1016/0010-4655(75)90055-7


% Definitions and initialization
MN1 = M + N + 1; MN2 = MN1 + 1;
M1 = M + 1; M2 = M + 2; M3 = M + 3;
N1 = N + 1; N2 = N + 2; N3 = N + 3;

% Initialize matrices A, B, and C
A = zeros(M1);
B = zeros(N1); B(1, 1) = 1;
CC = zeros(20, 20);
RH = zeros(20, 1);
NEQ = 0;

if (M < 0)
    error('M must be greater than or equal to 0.');
end

if (N <= 0)
    error('N must be greater than 0.');
end

%% PART 1 CALCULATION OF DENOMINATOR COEFFICIENTS (B-MATRIX)
% Definitions for loops
for IP = 0:N - 1 % last prong computed separately
    IP1 = IP + 1; IP2 = IP + 2;
    LP = N - IP; LP1 = LP + 1;
    LTP = M + N + 1 - IP; LTP1 = LTP + 1; LTP2 = LTP + 2;
    I = 0;
    IPS = M2;
    IPF = LTP;

    if (IP > M)
        IPS = IP2;
        IPF = N1;
    end

    % Label 2
    for ID = IPS:IPF % DO 9
        I = I + 1;

        %% Form RHS
        if IP == 0 % a) ZEROTH PRONG
            % a) 1. RH-VECTOR
            RH(I) = -C(1, ID);
            RH(I + N) = -C(ID, 1);
            %
            J = 0;
            % GO TO 7
        else % IP > 0,  b) INTERNAL PRONGS
            % Label 3
            % b) 1. RH-VECTOR BAR LAST ENTRY
            X1 = 0.0;
            X2 = 0.0;
            LT = ID;
            if ID > N1
                LT = N1;
            end
            for IS = 1:IP % DO 5
                for IT = 1:LT % DO 4
                    IPS2 = IP - IS + 2;
                    IDT1 = ID - IT + 1;
                    X1 = X1 + B(IS, IT) * C(IPS2, IDT1);
                    X2 = X2 + B(IT, IS) * C(IDT1, IPS2);
                end % Label 4
                IDS1 = ID - IS + 1;
                X1 = X1 + B(IP1, IS) * C(1, IDS1);
                X2 = X2 + B(IS, IP1) * C(IDS1, 1);
            end % Label 5
            RH(I) = -X1;
            ILP = I + LP;
            RH(ILP) = -X2;

            % b) 2. FIRST ROW OF CC-MATRIX
            % Label 6
            J = 1;
            IDP = ID - IP;
            CC(I, 1) = C(1, IDP);
            CC(ILP, 1) = C(IDP, 1);
        end

        % a) 2. TOTAL CC - MATRIX OF ZEROTH PRONG
        % b) 2. CC - MATRIX BAR FIRST & LAST ROWS OF INTERNAL PRONGS
        % Label 7
        for IT = IP2:N1 % DO 8
            J = J + 1;
            JLP = J + LP;
            ILP = I + LP;
            IDT1 = ID - IT + 1;
            CC(I, J) = 0.0;
            CC(I, JLP) = 0.0;
            CC(ILP, JLP) = 0.0;
            CC(ILP, J) = 0.0;

            if ID >= IT
                CC(I, J) = C(1, IDT1);
                CC(ILP, JLP) = C(IDT1, 1);
            end
        end % Label 8
    end % Label 9

    if IP > 0 % Else GO TO 17
        % b) 1. RH-VECTOR LAST ENTRY & b) 2. CC-MATRIX LAST ROW
        I = 2 * LP + 1;
        X = 0.0;

        if IP <= M % Else GO TO 13
            % b) 1. RH-VECTOR LAST ENTRY FOR TYPE (1) PRONGS
            for IS = 1:IP % DO 11
                for IT = 1:N1 % DO 10
                    IPS2 = IP - IS + 2;
                    LPT2 = LTP - IT + 2;
                    X = X + B(IS, IT) * C(IPS2, LPT2) + B(IT, IS) * C(LPT2, IPS2);
                end % Label 10
                LPS2 = LTP - IS + 2;
                X = X + B(IP1, IS) * C(1, LPS2) + B(IS, IP1) * C(LPS2, 1);
            end % Label 11
            RH(I) = -X;
            % b) 2. CC-MATRIX  LAST ROW FOR TYPE (1) PRONGS
            LT2P = LTP1 - IP;
            CC(I, 1) = C(1, LT2P) + C(LT2P, 1);
            for J = 2:LP1 % DO 12
                IT = J + IP;
                LPT2 = LTP - IT + 2;
                JLP = J + LP;
                CC(I, J) = C(1, LPT2);
                CC(I, JLP) = C(LPT2, 1);
            end % Label 12
            % GO TO 17
        else % GO TO 13
            % b) 1. RH-VECTOR LAST ENTRY FOR TYPE (2) PRONGS
            % Label 13
            for IS = 1:IP % DO 15
                for IT = 1:IP % DO 14
                    IPS2 = IP - IS + 2;
                    IPT2 = IP - IT + 2;
                    X = X + B(IS, IT) * C(IPS2, IPT2);
                end % Label 14
                X = X + B(IP1, IS) * C(1, IPS2) + B(IS, IP1) * C(IPS2, 1);
            end % Label 15
            RH(I) = -X;
            % b) 2. CC-MATRIX  LAST ROW FOR TYPE (2) PRONGS
            CC(I, 1) = C(1, 1);
            for J = 2:I % DO 16
                CC(I, J) = 0.0;
            end % Label 16
        end
    end

    %% Solve CC y = RH(Ax=b)
    % Label 17
    % INVERSION OF CC-MATRIX
    NEQ = 2 * LP + 1;
    if IP == 0
        NEQ = 2 * N;
    end

    % Call SOLVE function
    y = acrs_solve(CC, RH, NEQ);

    % Process results and reshape as B matrix
    % SOLUTION TO PRONG COEFFICIENTS (B-MATRIX)
    B(1, 1) = 1.0;
    if IP == 0
        for I = 1:N % DO 20
            B(1, I + 1) = y(I);
            B(I + 1, 1) = y(N + I);
        end % Label 20
    else
        B(IP + 1, IP + 1) = y(1);
        for I = 2:LP1 % DO 18
            IPI = IP + I;
            ILP = I + LP;
            B(IP1, IPI) = y(I);
            B(IPI, IP1) = y(ILP);
        end % Label 18
    end
    % INCREASE PRONG NUMBER AND REPEAT,UNLESS FINAL PRONG
end

% c) FINAL PRONG
% Label 22
X = 0.0;
if M >= N
    for IS = 1:N % DO 24
        for IT = 1:N % DO 23
            MS3 = M3 - IS;
            NT2 = N2 - IT;
            X = X + B(IS, IT) * C(MS3, NT2) + B(IT, IS) * C(NT2, MS3);
        end
        NS2 = N2 - IS;
        M2N = M2 - N;
        X = X + B(IS, N1) * (C(NS2, M2N) + C(MS3, 1)) + B(N1, IS) * (C(M2N, NS2) + C(1, MS3));
    end % Label 24
    Z = C(M2N, 1) + C(1, M2N);
    if abs(Z) < 1e-10
        % GO TO 31
        error('Failure on last prong.');
    end
    B(N1, N1) = -X / Z;
    % GO TO 28
else % M < N
    % Label 25
    for IS = 1:N % DO 27
        for IT = 1:N % DO 26
            NS2 = N2 - IS;
            NT2 = N2 - IT;
            X = X + B(IS, IT) * C(NS2, NT2);
        end % Label 26
        X = X + B(N1, IS) * C(1, NS2) + B(IS, N1) * C(NS2, 1);
    end % Label 27
    B(N1, N1) = -X;
end

% Label 28
%% PART 2 CALCULATION OF NUMERATOR COEFFICIENTS (A-MATRIX)
for IG = 1:M1 % DO 30
    for ID = 1:M1 % DO 30
        X = 0.0;
        LS = IG;
        if IG > N1
            LS = N1;
        end
        for IS = 1:LS % DO 29
            LT = ID;
            if ID > N1
                LT = N1;
            end
            for IT = 1:LT % DO 29
                IGS1 = IG - IS + 1;
                IDT1 = ID - IT + 1;
                X = X + B(IS, IT) * C(IGS1, IDT1);
            end % Label 29
        end % Label 29
        A(IG, ID) = X;
    end % Label 30
end % Label 30

% END function SOD()
end
