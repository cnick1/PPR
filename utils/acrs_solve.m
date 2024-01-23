function B = acrs_solve(A, B, N)
% Subroutine to solve SOD equations
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

% Initialize matrices and variables
M = 1;
DET = 1;
% Label 1
IPVOT = zeros(N, 1);
INDEX = zeros(N, 2);
PIVOT = zeros(N, 1);

% Search for pivot element
for I = 1:N % DO 24
    T = 0;
    for J = 1:N % DO 6
        if IPVOT(J) ~= 1
            for K = 1:N % DO 5
                if IPVOT(K) ~= 1
                    if abs(T) < abs(A(J, K))
                        IROW = J;
                        ICOL = K;
                        T = A(J, K);
                    end
                end
            end % Label 5
        end
    end % Label 6

    IPVOT(ICOL) = IPVOT(ICOL) + 1;

    % Put pivot element on diagonal
    if IROW ~= ICOL
        % Label 7
        DET = -DET;
        for L = 1:N
            temp = A(IROW, L);
            A(IROW, L) = A(ICOL, L);
            A(ICOL, L) = temp;
        end

        % if M > 0
        for L = 1:M
            temp = B(IROW, L);
            B(IROW, L) = B(ICOL, L);
            B(ICOL, L) = temp;
        end
        % end
    end

    % Label 11

    INDEX(I, 1) = IROW;
    INDEX(I, 2) = ICOL;
    PIVOT(I) = A(ICOL, ICOL);
    IM1 = I - 1;

    % Check for pivot too small
    if I ~= 1 && abs(PIVOT(I)) < exp(log(abs(DET)) / IM1) * I * 1e-10
        error(['PIVOT TOO SMALL, EXIT SINCE RANK=', num2str(IM1)]);
    end

    DET = DET * PIVOT(I);

    % Divide pivot row by pivot element
    A(ICOL, ICOL) = 1.0;
    for L = 1:N
        A(ICOL, L) = A(ICOL, L) / PIVOT(I);
    end

    % if M > 0
    for L = 1:M
        B(ICOL, L) = B(ICOL, L) / PIVOT(I);
    end
    % end

    % Reduce non-pivot rows
    for LI = 1:N
        if LI ~= ICOL
            T = A(LI, ICOL);
            A(LI, ICOL) = 0.0;
            for L = 1:N
                A(LI, L) = A(LI, L) - A(ICOL, L) * T;
            end

            % if M > 0
            for L = 1:M
                B(LI, L) = B(LI, L) - B(ICOL, L) * T;
            end
            % end
        end
    end
end % Label 24

end % End function acrs_solve
