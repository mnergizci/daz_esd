import os
import argparse
from tqdm import tqdm
import sys
def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Run daz_esd_01_prepare_input.py for each frame in the frame list.')
    parser.add_argument('frames_list', help='Path to the frame list file.')
    parser.add_argument('--add_eu', action='store_true', help='Flag to add EU')
    parser.add_argument('--s1ab', action='store_true', help='s1ab offset in step5')
    parser.add_argument('--roll_assist', action='store_true', help='rolling assistant for step5')
    
    ##start-stop section:
    parser.add_argument('--skip_step1', action='store_true', help='Skip step 1')
    parser.add_argument('--skip_step2', action='store_true', help='Skip step 2')
    parser.add_argument('--skip_step3', action='store_true', help='Skip step 3')
    parser.add_argument('--skip_step4', action='store_true', help='Skip step 4')
    parser.add_argument('--skip_step5', action='store_true', help='Skip step 5')
    
    args = parser.parse_args()

    # Get the batch directory from the environment variable or use './' as default
    batchdir = os.environ.get('BATCH_CACHE_DIR', './')
    datadir = os.path.join(batchdir, 'daz_esd')

    # Read the frame list file
    with open(args.frames_list, 'r') as file:
        lines = file.readlines()

    # First Round: Iterate over each line in the frame list with a progress bar
    if not args.skip_step1:
        for i, line in enumerate(tqdm(lines, desc="Processing frames", unit="frame")):
            # Strip any leading/trailing whitespace (including newlines)
            frame = line.strip()
            frame_step1_file = os.path.join(datadir, frame + '.step1.csv')
    
            if os.path.exists(frame_step1_file):
                print(f'{frame} already exists, skipping this frame')
            else:
                print(f'Frame {i+1}/{len(lines)}: {frame} processing!')
                print('step1!')
                # Define the command to run the Python script with the frame as an argument
                command = f"python daz_esd_01_prepare_input.py {frame} --orbdiff_fix --daz_wrt"
    
                # Execute the command
                os.system(command)

        print("First round of processing complete.")

    # Second Round: Check and run the second step for each frame. SET correction
    if not args.skip_step2:
        for line in lines:
            frame = line.strip()
            frame_step1_file = os.path.join(frame + '.step1.csv')
            frame_earthtides_file = os.path.join(frame + '.earthtides.csv')
            frame_step2_file = os.path.join(frame + '.step2.csv')
    
            if not os.path.exists(os.path.join(datadir,frame_step2_file)):
                print(f'{frame} step2 processing!')
                # Define the command to run the second Python script with the necessary arguments
                command = f"python daz_esd_02_extract_SET.py {frame_step1_file} {frame_earthtides_file} {frame_step2_file}"
                # Execute the command
                os.system(command)
            elif not os.path.exists(os.path.join(datadir,frame_earthtides_file)):
                # Handle the case if the .earthtides.csv file is missing (if needed)
                # For now, we just print a message
                print(f'Warning: Earthtides file for {frame} is missing!')
    
        print("Second round of processing complete.")

    
    # Third Round: Check and run the third step for each frame. Iono Extraction
    if not args.skip_step3:
        for line in lines:
            frame = line.strip()
            frame_step2_file = os.path.join(frame + '.step2.csv')
            frame_step3_file = os.path.join(frame + '.step3.csv')   ##it will used for slope calculationg in step5
    
            if not os.path.exists(os.path.join(datadir,frame_step3_file)):
                print(f'{frame} step3 processing!')
                # Define the command to run the third Python script with the necessary arguments
                command = f"python daz_esd_03_extract_iono.py {frame_step2_file} {frame_step3_file}"
                os.system(command)
                print("Third round of processing complete.")
            else:
                print('The Ionospheric correction already applied. Skip next step..')
                # Execute the command
    


    # Fourth Round: Check and run the fourth step for each frame. PPM extraction
    if not args.skip_step4:
        for line in lines:
            frame = line.strip()
            frame_step3_frame_file = os.path.join(frame + '.frame.step3.csv')
            frame_step4_frame_file = os.path.join(frame + '.frame.step4.csv')
            nc_file = os.path.join(datadir, 'vel_gps_kreemer.nc')
    
            if not os.path.exists(os.path.join(datadir,frame_step4_frame_file)):
                if args.add_eu:
                    print('add_eu active')
                    command = f'python daz_esd_04_extract_PMM.py {frame_step3_frame_file} {frame_step4_frame_file} {nc_file} --add_eu'
                else:
                    command = f'python daz_esd_04_extract_PMM.py {frame_step3_frame_file} {frame_step4_frame_file} {nc_file}'
                    print('no add_eu')
                # Execute the command
                os.system(command)
    
        print("Fourth round of processing complete.")


    
    # Fifth Round: Calculate slope and intercept using huber regression. 
    if not args.skip_step5:
        for line in lines:
            frame = line.strip()
            frame_step3_file=os.path.join(frame +'.step3.csv')
            frame_step4_frame_file=os.path.join(frame +'.frame.step4.csv')
            frame_step5_file=os.path.join(frame + '.step5.csv')
            frame_step5_frame_file=os.path.join(frame + '.frame.step5.csv')
    
            if not os.path.exists(os.path.join(datadir,frame_step5_file)) or not os.path.exists(os.path.join(datadir,frame_step5_frame_file)):
                if args.s1ab and args.roll_assist:
                    command= command= f'python daz_esd_05_calculate_slop.py {frame_step3_file} {frame_step4_frame_file} {frame_step5_file} {frame_step5_frame_file} --s1ab --roll_assist'
                elif args.s1ab:
                    print('s1ab offset active')
                    command= f'python daz_esd_05_calculate_slop.py {frame_step3_file} {frame_step4_frame_file} {frame_step5_file} {frame_step5_frame_file} --s1ab'
                elif args.roll_assist:
                    command= f'python daz_esd_05_calculate_slop.py {frame_step3_file} {frame_step4_frame_file} {frame_step5_file} {frame_step5_frame_file} --roll_assist'
                else:
                    command= f'python daz_esd_05_calculate_slop.py {frame_step3_file} {frame_step4_frame_file} {frame_step5_file} {frame_step5_frame_file}'
            
                # Execute the command
                os.system(command)
        print("Fifth round of processing complete.")

if __name__ == '__main__':
    main()
