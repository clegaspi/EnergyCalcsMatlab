classdef Tchain < handle
   % Simple site model of excitons, including field effects
   % Each site has a HOMO and LUMO with a gap Eo
   % adjacent HOMO's are coupled with betaH
   % adjacent LUMO's are coupled with betaE
   % electron and hole are coupled with a hubbard potential
   properties
      nsites     % number of sites on the chain
      nbasis     % number of s-ci basis functions
      e          % (1,nbasis) elec location for ibasis
      h          % (1,nbasis) hole location for ibasis
      periodic   % use periodic boundary conditions

      Eo         % HOMO LUMO gap on site
      aSite      % Length of a site in angstroms
      betaE      % coupling between adjacent LUMOs
      betaH      % coupling between adjacent HOMOs
      U          % e-e repulsion on a site
      eps        % dielectric constant
      
      dist         % distance between two sites (with or without periodicity)
      betaMatrixE  % multiplied by betaE and added to Hamiltonain
      betaMatrixH  % multiplied by betaH and added to Hamiltonain
      EHDiff       % elec-hole distance in site units
   end % properties
   methods
      function res = Tchain(nsites,periodic)
         res.nsites = nsites;
         res.periodic = periodic;
         res.createBasis;
         [res.betaMatrixE res.betaMatrixH] = res.calcBetaMatrices;
         res.EHDiff = res.calcEHDiff;
         res.dist = 1:nsites;
         if (periodic)
            res.dist = min(res.dist,nsites-res.dist);
         end             
      end
      function createBasis(obj)
         obj.nbasis = obj.nsites^2;
         ic = 0;
         for ih = 1:obj.nsites
            for ie = 1:obj.nsites
               ic = ic + 1;
               obj.h(ic) = ih;
               obj.e(ic) = ie;
            end
         end
      end
      function res = sciMatrix(obj, field)
         res = obj.Eo .* diag(ones(1,obj.nbasis)) + ...
            obj.betaE .* obj.betaMatrixE + ...
            obj.betaH .* obj.betaMatrixH + ...
            diag(-1.0 * obj.hubbard(obj.EHDiff .* obj.aSite)) + ...
            diag(-1.0 *field .* (obj.e-obj.h) .* obj.aSite);
      end
      function res = transVector(obj)
         res = obj.EHDiff == 0;
      end
      function [emat hmat] = calcBetaMatrices(obj)
         emat = zeros(obj.nbasis,obj.nbasis);
         hmat = zeros(obj.nbasis,obj.nbasis);
         for ib = 1:obj.nbasis
            for jb = ib:obj.nbasis
               de = abs(obj.e(ib)-obj.e(jb));
               dh = abs(obj.h(ib)-obj.h(jb));
               if (obj.periodic)
                  de = min(de,obj.nsites-de);
                  dh = min(dh,obj.nsites-dh);
               end
               if (((de==1) && (dh==0)))
                  emat(ib,jb) = 1;
                  emat(jb,ib) = 1;
               elseif((de==0) && (dh==1))
                  hmat(ib,jb) = -1;
                  hmat(jb,ib) = -1;
               end
            end
         end
      end
      function res = calcEHDiff(obj)
         res = abs(obj.e - obj.h);
         if (obj.periodic)
            res = min(res,obj.nsites-res);
         end
      end
      function res = hubbard(obj,r)
         res = 14.397./sqrt( (14.397./obj.U).*(14.397./obj.U).*ones(size(r)) ...
            + (obj.eps * r).^2);
      end
   end % methods
end %
