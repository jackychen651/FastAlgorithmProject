classdef LippmannSchwinger

    properties (SetAccess = public)
        GFFT
        nu
        ne
        me
        n
        m
        omega
        quadRule
    end
     
    methods
        function LS = LippmannSchwinger(x,y,omega,nu,a)
            n = length(x);
            m = length(y);
            X = repmat(x', 1, m);
            Y = repmat(y, n,1);
            Lp = 4*a ; 
            L  = 2*a;
            kx = (-(2*n):1:(2*n-1));
            ky = (-(2*m):1:(2*m-1));
            KX = (2*pi/Lp)*repmat(kx', 1, 4*m);
            KY = (2*pi/Lp)*repmat(ky, 4*n,1);

            S = sqrt(KX.^2 + KY.^2);

            G2D = @(L,k,s)(1 + ...
                           (1i*pi/2*L*besselh(0,1,L*k)).*(s.*besselj(1,L*s)) - ...
                           (1i*pi/2*L*k*besselh(1,1,L*k)).*besselj(0,L*s)...
                           )./(s.^2 - k^2);

            LS.GFFT = G2D(L, omega, S);

            LS.n = n;
            LS.m = m;
            LS.ne = 4*n;
            LS.me = 4*m;
            LS.omega = omega;
            
            nu_vect = nu(X,Y);
            LS.nu = nu_vect(:);

            if ~isempty(isnan(LS.nu))
                % fprintf("some values of the window may be Nan")
                LS.nu(isnan(LS.nu)) = 0.0;
            end
        end
        
        function Gu = apply_Green(LS, u)
            
            BExt = zeros(LS.ne,LS.me);
            BExt(1:LS.n,1:LS.m)= reshape(u,LS.n,LS.m);
            BFft = fftshift(fft2(BExt));
            BFft = LS.GFFT.*BFft;
            BExt = ifft2(ifftshift(BFft));
            Gu = BExt(1:LS.n, 1:LS.m);
        end
        
    end
end