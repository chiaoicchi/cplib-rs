{
  description = "competitive programming library";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    rust-overlay = {
      url = "github:oxalica/rust-overlay";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    cpbin-rs = {
      url = "github:chiaoicchi/cpbin-rs";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs =
    {
      self,
      nixpkgs,
      rust-overlay,
      cpbin-rs,
    }:
    let
      system = "x86_64-linux";
      pkgs = import nixpkgs {
        inherit system;
        overlays = [ rust-overlay.overlays.default ];
      };
      toolchain = pkgs.rust-bin.fromRustupToolchainFile ./rust-toolchain.toml;
      ckTool = pkgs.writeShellApplication {
        name = "ck";
        runtimeInputs = [
          pkgs.git
          pkgs.coreutils
          pkgs.findutils
          pkgs.diffutils
          pkgs.jq
          pkgs.python3
          pkgs.gcc
          toolchain
        ];
        text = builtins.readFile ./tools/ck.sh;
      };
      bdTool = pkgs.writeShellApplication {
        name = "bd";
        runtimeInputs = [
          pkgs.git
          pkgs.coreutils
          pkgs.gnused
          pkgs.rustfmt
          cpbin-rs.packages.${system}.bundle-rs
          pkgs.wl-clipboard
        ];
        text = builtins.readFile ./tools/bd.sh;
      };
    in
    {
      packages.${system} = {
        ck = ckTool;
        bd = bdTool;
      };
      devShells.${system}.default = pkgs.mkShell {
        packages = [
          (pkgs.rust-bin.fromRustupToolchainFile ./rust-toolchain.toml)
          ckTool
          bdTool
        ];
        shellHook = ''
          echo "cplib environment"
          echo "  rust: $(rustc --version)"
        '';
      };
    };
}
