# Introduction
This is a SEGMENT-specific fork of scE2G-pipeline, which consists of 2 (nested) submodules which are in turn also forks of the original repositories (ENCODE-rE2G and ABC).

After making changes to either of the submodules, commit and push changes in an inside-out fashion, so the pointers in the outer repositories are updated.
This can be implemented more easy using the following helper-function, which can be added to your `.bashrc` or `.zshrc` file.

```bash
rgcmp () {

  orig_dir=$(pwd)
  commit_message="$1"

  # Commit changes in the current repository if there are any
  if [ -n "$(git status --porcelain)" ]; then
      git commit -m "$commit_message"
      git push
  else
      echo "------------No changes to commit in $(pwd)------------"
  fi

  # While the current repository is a submodule (i.e. it has a superproject)
  while parent=$(git rev-parse --show-superproject-working-tree 2>/dev/null) && [ -n "$parent" ]; do
      # Use realpath to get the relative path from the parent's working tree to our current directory.
      current_dir=$(pwd)
      relative_path=$(realpath --relative-to="$parent" "$current_dir")

      # Change directory to the parent repository
      builtin cd "$parent" || exit 1

      # Stage the submodule pointer change (if there is one)
      if [ -n "$(git status --porcelain "$relative_path")" ]; then
          echo "------------Updating pointer in $relative_path------------"
          git add "$relative_path"
          git commit -m "Update submodule '$relative_path' pointer"
          git push
      else
          echo "------------No pointer update to commit for '$relative_path' in parent repository------------"
      fi
  done

  builtin cd $orig_dir

}
```

# Notable changes from original include:
- Revised determine_mem function with scaler argument which originally set to 4, but can now be set per rule, specifically used in compute_kendall.smk
- IGV typo fix (issue supplied, will likely end up in upstream)
- Removed unncessary output_files in workflow/Snakefile
