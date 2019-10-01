for name in 2019-*.txt 
do
    n=${name%.txt}    # name without filename extension
    first=${n%% -0.*}  # first part of name
    last=${n#* -0.}    # last part of name
    new="$0.last-$first.txt"  # recombined new name

    printf 'Would move "%s" to "%s"\n' "$name" "$new"
    # mv -- "$name" "$new"
done

# run the following script from a directory
for d in ./*/; do (cd "$d" && for name in 2019-*.txt; do n=${name%%.txt}; first=${n%%-0.*}; last=${n#2019-*-0}; new="0$last-$first"; echo "$name to $new"; mv -i "$name" "$new"; done; ); done;
