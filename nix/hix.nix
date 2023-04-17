{pkgs, ...}: {
  name = "impute";
  compiler-nix-name = "ghc944"; # Version of GHC to use

  # crossPlatforms = p: pkgs.lib.optionals pkgs.stdenv.hostPlatform.isx86_64 ([
  #   p.mingwW64
  # ] ++ pkgs.lib.optionals pkgs.stdenv.hostPlatform.isLinux [
  #   p.musl64
  # ]);

  # Tools to include in the development shell
  shell.tools.cabal = "latest";
  #shell.tools.haskell-language-server = "latest";
}
