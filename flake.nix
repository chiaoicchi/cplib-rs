{
  description = "competitive programming library";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    rust-overlay = {
      url = "github:oxalica/rust-overlay";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    bundle-rs = {
      url = "github:chiaoicchi/bundle-rs";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.rust-overlay.follows = "rust-overlay";
    };
  };

  outputs =
    {
      self,
      nixpkgs,
      rust-overlay,
      bundle-rs,
    }:
    let
      system = "x86_64-linux";
      pkgs = import nixpkgs {
        inherit system;
        overlays = [ rust-overlay.overlays.default ];
      };
      toolchain = pkgs.rust-bin.fromRustupToolchainFile ./rust-toolchain.toml;
      bundleRs = bundle-rs.packages.${system}.bundle-rs;
      bdTool = bundle-rs.packages.${system}.bd;
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
          bundleRs
        ];
        text = builtins.readFile ./tools/ck.sh;
      };
    in
    {
      packages.${system} = {
        ck = ckTool;
      };
      devShells.${system}.default = pkgs.mkShell {
        packages = [
          toolchain
          ckTool
          bdTool
          bundleRs
        ];
        shellHook = ''
          export CPLIB_ROOT="$(git rev-parse --show-toplevel)"
          echo "cplib environment"
          echo "  rust: $(rustc --version)"
          echo "  bundle-rs: $(command -v bundle-rs)"
        '';
      };
    };
}
